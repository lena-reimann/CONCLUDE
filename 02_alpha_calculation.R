#######################
## Alpha calculation ##
#######################
# by Lena Reimann
# Nov 1, 2024

# Goal: calculate A for each calibration period to correct for local processes that make the pop results deviate from the results solely based on the gravity approach

rm(list=ls())

lib = "pathtolibraries"

library(sp, lib.loc = lib)
library(raster, lib.loc = lib)
library(geosphere, lib.loc = lib)
library(scales, lib.loc = lib) 
library(Matrix, lib.loc = lib)

#for parallel processing
library(parallel, lib.loc = lib) 
library(foreach, lib.loc = lib) 
library(doParallel, lib.loc = lib) 
library(igraph, lib.loc = lib) 
library(doSNOW, lib.loc = lib)


#-------------------------------------------#
# Step 1 - Define basic information of task #
#-------------------------------------------#
task         <- "alpha_calc"
reg          <- "N"              # change accordingly
iso           <- c("SVN")   # change accordingly
id            <- c(705) 
p             <- "P1"             # change accordingly
gw            <- ifelse(reg == "SE", 10, 20)

settl           <- c("r", "u")

cores <- 4

# paths
path <- "path_to_data"
setwd(paste(path, sep = "/")) 
getwd()

path_m = "path_to_distance_matrix"

#-----------------------------------------------------#
# Step 2 - define function for alpha calc and scaling #
#-----------------------------------------------------#

# (a) calculate alphas for cst versus inl locations (run for urb and rur individually)
alpha_calc <- function(beta_cst, beta_inl, vec_alpha, cells_dist, cst_mask, inl_mask, base_year_totpop,
                       cntr_li, cntr_ll, chg_cst_base_to_next, chg_inl_base_to_next, pop_next, pop_base) {
  
  # (1) calculate pje and potential for cst versus inl
  # calculate pje cst
  exp_cells_dist <- cells_dist
  exp_cells_dist@x <- exp(-beta_cst * exp_cells_dist@x)
  pje_cst <- as.vector(exp_cells_dist %*% base_year_totpop) 
  
  # calculate pje inl
  exp_cells_dist <- cells_dist
  exp_cells_dist@x <- exp(-beta_inl * exp_cells_dist@x)
  pje_inl <- as.vector(exp_cells_dist %*% base_year_totpop) 
  
  # calculate potentials
  pot_cst <- cntr_li * (pje_cst + ((base_year_totpop * cst_mask) * vec_alpha))
  pot_inl <- cntr_li * (pje_inl + ((base_year_totpop * inl_mask) * vec_alpha))
  
  # (2) calculate alphas
  # calculate v_hat
  v_hat_cst <- (((pop_next * cst_mask) - (pop_base * cst_mask)) / chg_cst_base_to_next) * sum(pot_cst)
  v_hat_inl <- (((pop_next * inl_mask) - (pop_base * inl_mask)) / chg_inl_base_to_next) * sum(pot_inl)
  
  # alphas
  alpha_cst <- (v_hat_cst / cntr_li - pje_cst) / (base_year_totpop * cst_mask)
  alpha_inl <- (v_hat_inl / cntr_li - pje_inl) / (base_year_totpop * inl_mask)
  
  # postprocess results
  alpha_cst[is.na(alpha_cst)] <- 0
  alpha_inl[is.na(alpha_inl)] <- 0
  
  alpha_cst[!is.finite(alpha_cst)] <- 0
  alpha_inl[!is.finite(alpha_inl)] <- 0  
  
  # make df
  alpha <- data.frame(cntr_ll, (alpha_cst + alpha_inl))
  colnames(alpha) [3] <- "alpha"
  
  return(alpha)
}

# (b) scale alphas
scale_alpha <- function(alpha, quantiles) { 
  
  # 1. calculate quantiles
  qu_neg <- quantile(alpha[which(alpha[,"alpha"]<0),"alpha"], probs = c(.05, .10, .25, .33, .50, .66, .75, .83, .90, .95))
  qu_pos <- quantile(alpha[which(alpha[,"alpha"]>0),"alpha"], probs = c(.05, .10, .25, .33, .50, .66, .75, .83, .90, .95))
  
  if (quantiles == "a_q50") {
    # remove upper and lower 25%
    alpha[which(alpha[,"alpha"]> qu_pos["75%"]),"alpha"] <- qu_pos["75%"] 
    alpha[which(alpha[,"alpha"]< qu_neg["25%"]),"alpha"] <- qu_neg["25%"] 
  } else {
    # remove upper and lower 10%
    alpha[which(alpha[,"alpha"]> qu_pos["90%"]),"alpha"] <- qu_pos["90%"] 
    alpha[which(alpha[,"alpha"]< qu_neg["10%"]),"alpha"] <- qu_neg["10%"] 
  }
  
  # 2. scale alphas from 0-100
  max <- max(alpha[,"alpha"])
  min <- min(alpha[,"alpha"])
  # establish "max" for scaling factor
  max <- ifelse(max > abs(min), max, abs(min))
  # calc scaling factor
  factor <- 100 / max
  # scale alphas (preserve skewness of the function)
  alpha$alpha_sc <- alpha$alpha * factor
  
  return(alpha)
  
  }


#------------------------#
# Step 3 - Load the data #
#------------------------#

# load country ids based on Gridded Population of the World (GPW)
gpw_ids <- read.csv("cntr_ids.csv")
isos <- gpw_ids$iso
ids <- gpw_ids$gpw_id

# pop data
t = (1990, 2000, 2010)
pop_rur <- stack(paste0("rur_pop", t, ".tif"))
pop_urb <- stack(paste0("urb_pop", t, ".tif"))  
pop_tot <- stack(paste0("tot_pop", t, ".tif"))

# spatial mask (for all countries)
mask00 <- stack(paste0(path, "/sp_mask_", ids, "_2000.tif"))

# defining coastal areas (= Low Elevation Coastal Zone (LECZ) + 20 km coastline buffer)
cst <-raster(paste(path, "cst.tif", sep = "/")) 

# betas
# see beta cali insights for more info on betas used; adjusted based on urban sprawl analysis
if (reg == "SE") {
  if (p == "P1") { 
    beta_cr <- 0.004875
    beta_cu <- 0.543511617
    beta_ir <- 0.005125
    beta_iu <- 0.729270887
    } else if (p == "P2") {
      beta_cr <- 0.184926178
      beta_cu <- 0.539798844
      beta_ir <- 0.204392091
      beta_iu <- 0.659754143
      } else {
        beta_cr <- 0.0840616795915027 # averaged for P12
        beta_cu <- 0.524418163450559
        beta_ir <- 0.0929102774432398
        beta_iu <- 0.709506927021345
  }
} else {
  if (p == "P1") {
    beta_cr <- 1.82 # calibrated for P1
    beta_cu <- 0.344505198
    beta_ir <- 2.18
    beta_iu <- 0.380768903
    } else if (p == "P2") {
      beta_cr <- 1.82  # calibrated for P2
      beta_cu <- 0.500517289
      beta_ir <- 2.18
      beta_iu <- 0.611743354
      } else {
        beta_cr <- 1.8 # averaged for P12 and adjusted based on urban sprawl analysis
        beta_cu <- 0.43641450168164
        beta_ir <- 2.2
        beta_iu <- 0.482352870279708
        }
}


#----------------------------#
# preprocessing + alpha calc #
#----------------------------#

# select data for time steps
if (p == "P1") {
  pop_rur <- pop_rur[[1:2]]
  pop_urb <- pop_urb[[1:2]]
  pop_tot <- pop_tot[[1:2]]
} else if (p == "P2") {
  pop_rur <- pop_rur[[2:3]]
  pop_urb <- pop_urb[[2:3]]
  pop_tot <- pop_tot[[2:3]]
} else {
  pop_rur <- pop_rur[[c(1,3)]]
  pop_urb <- pop_urb[[c(1,3)]]
  pop_tot <- pop_tot[[c(1,3)]]
}

# lists for calculated alphas as tifs
alpha_u_sc_tif <- list()
alpha_r_sc_tif <- list()

for (i in 1:length(id)) {
  
  ## (0) load distance matrix
  dist_m <- readMM(paste(paste(path_m, iso[i], sep = "/"), "mtx.gz", sep = "."))
  
  ## (1) select mask for cntr
  mask00_c <- mask00[[which(ids == id[i])]]  
  
  ## (2) preprocess pop data
  # crop to mask
  pop_rur_c <- crop(pop_rur, mask00_c)
  pop_urb_c <- crop(pop_urb, mask00_c)
  pop_tot_c <- crop(pop_tot, mask00_c)
  
  # set NAs to 0
  pop_rur_c[is.na(pop_rur_c)] <- 0
  pop_urb_c[is.na(pop_urb_c)] <- 0
  pop_tot_c[is.na(pop_tot_c)] <- 0
  
  # mask pop with mask
  pop_rur_c <- mask(pop_rur_c, mask00_c)
  pop_urb_c <- mask(pop_urb_c, mask00_c)
  pop_tot_c <- mask(pop_tot_c, mask00_c)
  
  pop_tot_c_mask <- pop_tot_c # mask for rasterizing later
  
  # convert to points
  pop_rur_c <- rasterToPoints(pop_rur_c)
  pop_urb_c <- rasterToPoints(pop_urb_c)
  pop_tot_c <- rasterToPoints(pop_tot_c)
  
  # define lat long and spatial mask values
  cntr_LL <- rasterToPoints(mask00_c)[,c(1,2)]
  cntr_Li <- rasterToPoints(mask00_c)[,3]
  
  ## (3) preprocess coastal mask
  cst_c <- crop(cst, mask00_c)
  cst_c <- mask(cst_c, mask00_c)
  
  inl_c <- 1-cst_c
  
  # convert to points
  cst_c <- rasterToPoints(cst_c)
  inl_c <- rasterToPoints(inl_c)

  ## (4) establish pop per region
  cr_b <- pop_rur_c[,3] * cst_c[,3]
  cr_e <- pop_rur_c[,4] * cst_c[,3]
  
  cu_b <- pop_urb_c[,3] * cst_c[,3]
  cu_e <- pop_urb_c[,4] * cst_c[,3]
  
  cst_b <- pop_tot_c[,3] * cst_c[,3]
  cst_e <- pop_tot_c[,4] * cst_c[,3]
    
  ir_b <- pop_rur_c[,3] * inl_c[,3]
  ir_e <- pop_rur_c[,4] * inl_c[,3]
  
  iu_b <- pop_urb_c[,3] * inl_c[,3]
  iu_e <- pop_urb_c[,4] * inl_c[,3]
  
  inl_b <- pop_tot_c[,3] * inl_c[,3]
  inl_e <- pop_tot_c[,4] * inl_c[,3]
  
  urb_e <- pop_urb_c[,4]
  urb_b <- pop_urb_c[,3]
  
  rur_e <- pop_rur_c[,4]
  rur_b <- pop_rur_c[,3]
  
  tot_e <- pop_tot_c[,4]
  tot_b <- pop_tot_c[,3]
  
  # calculate population change per region
  cu_chg <- sum(cu_e) - sum(cu_b)
  cr_chg <- sum(cr_e) - sum(cr_b)
  
  iu_chg <- sum(iu_e) - sum(iu_b)
  ir_chg <- sum(ir_e) - sum(ir_b)
  
  # check whether the sum of the rural and urban is equal to the total population
  #sum(cst_b) - (sum(cr_b) + sum(cu_b))
  #sum(inl_b) - (sum(ir_b) + sum(iu_b))
  
  #sum(cst_e) - (sum(cr_e) + sum(cu_e))
  #sum(inl_e) - (sum(ir_e) + sum(iu_e))
  
  ## run alpha_cal function ##
  alpha <- rep(1, length(cntr_Li))

  # rural alphas
  alpha_r <- alpha_calc(beta_cr,
                        beta_ir,
                        alpha,
                        dist_m,
                        cst_c[,3],
                        inl_c[,3],
                        tot_b,
                        cntr_Li,
                        cntr_LL,
                        cr_chg,
                        ir_chg,
                        rur_e,
                        rur_b
  )  
  # set rural alphas to 0 when urban pop at base year >0
  alpha_r[urb_b>0, 3] <- 0
  
  # urban alphas
  alpha_u <- alpha_calc(beta_cu,
                        beta_iu,
                        alpha,
                        dist_m,
                        cst_c[,3],
                        inl_c[,3],
                        tot_b,
                        cntr_Li,
                        cntr_LL,
                        cu_chg,
                        iu_chg,
                        urb_e,
                        urb_b
  ) 
  

  ## run scale_alpha function ##
  qu <- "a_q50"
  alpha_r <- scale_alpha(alpha_r, qu)
  alpha_u <- scale_alpha(alpha_u, qu)
  
  # write csvs
  write.csv(alpha_r, paste0(paste(paste(path, "chunks", iso[i], sep = "/"), "alpha_r", p, sep = "_"), ".csv"), 
            row.names = F)
  write.csv(alpha_u, paste0(paste(paste(path, "chunks", iso[i], sep = "/"), "alpha_u", p, sep = "_"), ".csv"), 
            row.names = F)
  
  
  #make tifs
  alpha_u_sc_tif[[i]] <- rasterize(cntr_LL, pop_tot_c_mask, alpha_u[,"alpha_sc"])
  alpha_r_sc_tif[[i]] <- rasterize(cntr_LL, pop_tot_c_mask, alpha_r[,"alpha_sc"])
  
} #end of i loop


#### -------------- run only if more countries used ------------------ ####

#### mosaic all alpha rasters ####
## make mosaic of scaled alphas for entire region --> needed for population projections
alpha_u_sc_tif$fun <- mean
alpha_r_sc_tif$fun <- mean

alpha_u_sc_tif$na.rm <- TRUE
alpha_r_sc_tif$na.rm <- TRUE

alpha_u_sc_tif <- do.call(mosaic, alpha_u_sc_tif)
alpha_r_sc_tif <- do.call(mosaic, alpha_r_sc_tif)

# write raster
writeRaster(alpha_u_sc_tif, paste0(paste(p, "alpha", "u", sep = "_"), ".tif"))
writeRaster(alpha_r_sc_tif, paste0(paste(p, "alpha", "r", sep = "_"), ".tif"))
