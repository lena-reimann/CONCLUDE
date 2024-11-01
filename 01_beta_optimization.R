###########################
## Run beta optimization ##
###########################
# by Lena Reimann
# Nov 1, 2024

# Goal: search the optimized value of beta for rural and urban population

rm(list=ls())

lib = "pathtolibraries"

library(raster, lib.loc = lib)
library(sf, lib.loc = lib)
library(tidyverse, lib.loc = lib)
library(geosphere, lib.loc = lib)
library(Matrix, lib.loc = lib)

library(parallel, lib.loc = lib) 
library(foreach, lib.loc = lib) 
library(doParallel, lib.loc = lib) 
library(igraph, lib.loc = lib) 
library(doSNOW, lib.loc = lib)

#-------------------------------------------#
# Step 1 - Define basic information of task #
#-------------------------------------------#
task            <- "beta_opt"
reg             <- "N"
iso             <- c("SVN")
id              <- c(705)
p               <- "P2"  # which calibration period? P1 = 1990-2000; P2 = 2000-2010; P12 = 1990-2010
mode            <- "tot" # differentiating cr, cu, ir, iu (= cst/rur) or urb, rur (= reg)?
gw              <- ifelse(reg == "N", 20, 10)

# needed for running the optimization function
cores           <- 4 #8 on cluster
iter            <- 100
alpha           <- 1

settl           <- c("r", "u")
tolerance       <- 1e-3

# beta ranges to explore during the optimization
range_lower    <- c(-2: 0)
range_upper    <- c( 1: 2)

lrange_lower   <- length(range_lower)
lrange_upper   <- length(range_upper)

# paths
path <- "path_to_data"
setwd(paste(path, sep = "/")) 
getwd()

path_m = "path_to_distance_matrix"

#---------------------------------------#
# Step 2 - define optimization function #
#---------------------------------------#

model_square_error_sparse <- function(beta, vec_alpha, cells_dist, base_year_totpop, cntr_li, change_base_to_next, pop_next, pop_base) {
  # Exponentiate non zero distances:
  exp_cells_dist <- cells_dist
  exp_cells_dist@x <- exp(-beta * exp_cells_dist@x)
  # Perform the same calculation with dense function (here no non_zero_dist_mask needed)
  base_year_potentials <- as.vector(exp_cells_dist %*% base_year_totpop) + base_year_totpop * vec_alpha
  base_year_potentials <- cntr_li * base_year_potentials
  
  # Calculate model's population according to the potentials
  if (change_base_to_next >= 0) {
    change <- base_year_potentials / sum(base_year_potentials) # CAREFUL: POSSIBLE DIVISION BY ZERO
    model_pop_next <- pop_base + change * change_base_to_next
  }
  else {
    pot <- base_year_potentials
    pop_itr <- pop_base
    pop_loss <- change_base_to_next

    iter <- 100
    for (j in seq_len(iter)) {
      r_decr             <- pot / sum(pot, na.rm = T)
      pop_itr            <- pop_itr + pop_loss * r_decr
      pop_loss           <- abs(sum(pop_itr[which(pop_itr < 0)]))
      
      pop_itr[which(pop_itr < 0)] <- 0 # pop for next iteration
      pot[pop_itr==0]         <- 0 # pot for next iteration
      
    
      if (abs(pop_loss) < 0.000001) {
        model_pop_next  <- pop_itr
        break
      }
    }
  }
  
  #Step 3. Difference btw Model & observed
  dif_m_obs  <- model_pop_next - pop_next
  sum_sq     <- sum(dif_m_obs^2)
  # step 4. return the result
  return(sum_sq)
}


#----------------------------------------#
# Step 3a - Load & preprocess data needed #
#----------------------------------------#

# define time steps
if (p == "P1") {
  t_pop1 <- c(1990, 2000)
  t_pop2 <-  c("90", "00")
} else if (p == "P2") {
  t_pop1 <- c(2000, 2010)
  t_pop2 <-  c("00", "10")
} else {
  t_pop1 <- c(1990, 2010)
  t_pop2 <-  c("90", "10")
}

# list for calibrated betas
optm_list <- list()

# loop through all countries to calibrate
for (j in 1:length(id)) {
  
  # load distance matrix
  dist_m <- readMM(paste(paste(path_m, iso[k], sep = "/"), "mtx.gz", sep = "."))
  
  # country shape to be used
  mask <- raster(paste0(paste(paste0(paste(paste(path, ssp, "sp_mask", sep = "/"), ssp, "id", sep = "_"), ids[i]), "2000", sep = "_"), ".tif") )

  # pop data
  #rural
  rur <- stack(paste0(path, "/pop_", t_pop1, "/pop", t_pop2, "_cntr25_rur/pop", t_pop2, "_R2018_rur_id", id[j], ".tif"))
  cr <- stack(paste0(path, "/pop_", t_pop1, "/pop", t_pop2, "_cntr25_cr/pop", t_pop2, "_R2018_cr_id", id[j], ".tif"))
  ir <- stack(paste0(path, "/pop_", t_pop1, "/pop", t_pop2, "_cntr25_ir/pop", t_pop2, "_R2018_ir_id", id[j], ".tif"))
  
  #urban 
  urb <- stack(paste0(path, "/pop_", t_pop1, "/pop", t_pop2, "_cntr25_urb/pop", t_pop2, "_R2018_urb_id", id[j], ".tif"))
  cu <- stack(paste0(path, "/pop_", t_pop1, "/pop", t_pop2, "_cntr25_cu/pop", t_pop2, "_R2018_cu_id", id[j], ".tif"))
  iu <- stack(paste0(path, "/pop_", t_pop1, "/pop", t_pop2, "_cntr25_iu/pop", t_pop2, "_R2018_iu_id", id[j], ".tif")) 
  
  #total
  tot <- stack(paste0(path, "/pop_", t_pop1, "/pop", t_pop2, "_cntr25/pop", t_pop2, "_R2018_cntr25_id", id[j], ".tif"))
  cst <- stack(paste0(path, "/pop_", t_pop1, "/pop", t_pop2, "_cntr25_cst/pop", t_pop2, "_R2018_cst_id", id[j], ".tif"))
  inl <- stack(paste0(path, "/pop_", t_pop1, "/pop", t_pop2, "_cntr25_inl/pop", t_pop2, "_R2018_inl_id", id[j], ".tif"))
  
  
  # replace NAs with 0s (otherwise not included in the points in the next step)
  #coastal
  cr[is.na(cr)] <- 0 
  cu[is.na(cu)] <- 0 
  #cst[is.na(cst)] <- 0 
   
  #inland
  ir[is.na(ir)] <- 0 
  iu[is.na(iu)] <- 0 
  #inl[is.na(inl)] <- 0 

  #total
  rur[is.na(rur)] <- 0 
  urb[is.na(urb)] <- 0 

  # define raster extent
  e <- extent(cr)
  
  # mask with gpw to reduce number of cells
  #mask <- crop(mask, e)
  
  #coastal
  cr <- mask(cr, mask)
  cu <- mask(cu, mask)
  #cst <- mask(cst, mask)

  #inland
  ir <- mask(ir, mask)
  iu <- mask(iu, mask)
  #inl <- mask(inl, mask)

  #total
  rur <- mask(rur, mask)
  urb <- mask(urb, mask)
  
  
  ## ------ V4 include values that are reclassified to urban in the rural calibration ------ ##
  cr[[2]][cr[[1]]>0 & cu[[2]]>0] <- cu[[2]][cr[[1]]>0 & cu[[2]]>0]
  ir[[2]][ir[[1]]>0 & iu[[2]]>0] <- iu[[2]][ir[[1]]>0 & iu[[2]]>0] 
   
  # convert the data to vectors
  #coastal
  cr <- rasterToPoints(cr)
  cu <- rasterToPoints(cu)
  
  #inland
  ir <- rasterToPoints(ir)
  iu <- rasterToPoints(iu)
  #inl <- rasterToPoints(inl)
 
  rur <- rasterToPoints(rur)
  urb <- rasterToPoints(urb)

  # convert mask to points
  mask <- rasterToPoints(mask)
  
  # take out the dev. potential
  cntr_Li <- mask[,3] 
  
  # take out the lat & lon
  cntr_LL <- mask[,c(1,2)]
  
  # make vectors of pop data
  cr_b <- cr[, 3] # first time step
  cr_e <- cr[, 4] # second time step
  
  cu_b <- cu[, 3]
  cu_e <- cu[, 4]

  ir_b <- ir[, 3]
  ir_e <- ir[, 4]
  
  iu_b <- iu[, 3]
  iu_e <- iu[, 4]
  
  rur_b <- rur[, 3]
  rur_e <- rur[, 4]
  
  urb_b <- urb[, 3]
  urb_e <- urb[, 4]
  
  tot_b <- rur_b + urb_b
  tot_e <- rur_e + urb_e
  
  # check whether the sum of the rural and urban is equal to the total population
  sum(rur_b) - (sum(ir_b) + sum(cr_b))
  sum(urb_b) - (sum(iu_b) + sum(cu_b))
  # ----
  sum(rur_e) - (sum(cr_e) + sum(ir_e))
  sum(urb_e) - (sum(cu_e) + sum(iu_e))
  # ----
  
  # population change
  totchang <- sum(tot_e) - sum(tot_b)
  rurchang <- sum(cr_e + ir_e) - sum(cr_b + ir_b)
  urbchang <- sum(cu_e + iu_e) - sum(cu_b + iu_b)
  
  # check sums
  totchang - (rurchang + urbchang)
  
  # caluclate vectors needed for optimization function
  rur_b <- cr_b + ir_b
  rur_e <- cr_e + ir_e 
  
  urb_b <- cu_b + iu_b
  urb_e <- cu_e + iu_e 
  
  
  #------------------------------------------------#
  # Step 3b - Run function for minimizing the error #
  #------------------------------------------------#
  num_cells <- length(tot_b)
  alphas <- c(rep(alpha, num_cells))
  
  
  # loop through betas in parallel
  cl <- makeCluster(cores)
  registerDoParallel(cl, cores = cores)
  
  t_start <- Sys.time()
  optm_tot <- 
    foreach(m = 1:length(settl), .combine = 'rbind') %:%
    foreach(i = 1:lrange_lower, .combine = 'rbind') %:%
    foreach(k = 1:lrange_upper, 
            .packages = c("raster", "igraph", "geosphere"),
            .combine = 'rbind') %dopar% {
              
              if (settl[m]=="u") { # urban
                optm <- optimize(model_square_error_sparse, c(range_lower[i], range_upper[k]), 
                                 vec_alpha = alphas,
                                 cells_dist = dist_m,
                                 base_year_totpop = tot_b,
                                 cntr_li = cntr_Li,
                                 change_base_to_next = urbchang,
                                 pop_next = urb_e,
                                 pop_base = urb_b,
                                 tol = tolerance)   
              } else { # rural
                optm <- optimize(model_square_error_sparse, c(range_lower[i], range_upper[k]), 
                                 vec_alpha = alphas,
                                 cells_dist = dist_m,
                                 base_year_totpop = tot_b,
                                 cntr_li = cntr_Li,
                                 change_base_to_next = rurchang,
                                 pop_next = rur_e,
                                 pop_base = rur_b,
                                 tol = tolerance)   
              }
              
              c(settl[m], range_lower[i], range_upper[k], optm[[1]], optm[[2]])
              
            } # end of foreach loops
  
  stopCluster(cl) 
  t_stop <- Sys.time()
  print(t_stop - t_start)
  
  optm_tot <- as.data.frame(optm_tot, row.names = F)
  colnames(optm_tot) <- c("settl", "lower", "upper", "betas", "err_sq")
  
  
  # write all results in opt list
  optm_list[[j]] <- optm_tot
  

} # end of j loop

## save workspace ##
save.image(paste0(path, task, "/", task, "_", p, ".RData"))



