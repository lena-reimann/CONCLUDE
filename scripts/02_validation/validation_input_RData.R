###########################
## Validation input data ##
###########################
# by Lena Reimann
# Jan 5, 2021


rm(list=ls())

#load packages
#library(rgdal)
library(raster)
#library(tmaptools)
library(geosphere)

#for parallel processing
#library("parallel")
#library("foreach")
#library("doParallel")
#library("igraph")
#library(doSNOW)

#cores <- 12

#--------------------------------------------------------------
# Section I  --  This Section Should be Modified Accordingly   
#--------------------------------------------------------------
# Make changes according to the specific country that we'd like to work on
#--------------------------------------------------------------
task          <- "validation"
ssp           <- "ssp3"     # change accordingly
b             <- "beta_v4_pos"


path          <- "D:/PhD/Fulbright/Demographic_model/model_extension" # change accordingly
setwd(paste(path, sep = "/"))       # may keep this line unchanged (set working directory) 
getwd()


# define time steps
t_pop1 <- c(1990, 2000, 2010, 2015)
t_pop2 <-  c("90", "00", "10", "15")
t_pop3 <-  c("pop90", "pop00", "pop10", "pop15")

#t <- seq(2010, 2100, by=10)  # projections years including 2010 (baseline)
#t1 <- seq(2000, 2100, by=10) # projections years including 2000 (baseline used for model calibration)
#t2 <- c(2000, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100) # for SSP1, leaving out the mask of 2010 because it has been produced based on SLR data already; the mask of 2000 is representative of the 2010 mask without SLR
 
#t_proj <- seq(2020,2100, by= 10) # maybe add fove year step later in case the process is fast (but: SLR variables also only produced in 10-year steps)


# define gpw ids
gpw_ids <- read.csv("boundary_data/gpw4_10/gpw_id_med_iso25.csv")
isos <- gpw_ids$iso
ids <- gpw_ids$gpw_id


# remove dza from list
#isos <- isos[(isos != "DZA")]
#ids <- ids[ids != 12]


# add established loss of attractiveness factors in ESL locations
#if (reg == "SE") {
#  cor_esl_u <- 1.30864407
#  cor_esl_r <- 0.50891803
#} else {
#  cor_esl_u <- 1.60596192
#  cor_esl_r <- 0.20403168
#}


#--------------------------------------------------------------
# Section II  --  get the data
#--------------------------------------------------------------
# (1) load mask data 
# for ssp3: mask of 2000

#entire med
mask00_med <- readAll(raster("spatial_mask/final/sp_mask_ssp3v5.tif"))

#per cntr
mask_cntr <- list()

# define t
t_sp_mask <- c(2000)

#load masks per country for each t
for (i in seq(ids)) {
  name <- paste0(paste(paste0(paste(paste0("spatial_mask/sp_mask_cntr25_", ssp, "/sp_mask"), ssp, "id", sep = "_"), ids[i]), 
                       t_sp_mask, sep = "_"), ".tif")  
  mask_cntr[[i]] <- readAll(stack(name))
}


# (2) load alphas
alpha00_u <- list()
alpha10_u <- list()

alpha00_r <- list()
alpha10_r <- list()

for (i in 1:length(ids)) { 
    alpha00_u[[i]] <- readAll(raster(paste0("calibration/slurm_NEW/chunks/results/region/", isos[i], "_P1_alpha_u_", b, ".tif")))
    alpha10_u[[i]] <- readAll(raster(paste0("calibration/slurm_NEW/chunks/results/region/", isos[i], "_P12_alpha_u_", b, ".tif")))
    
    alpha00_r[[i]] <- readAll(raster(paste0("calibration/slurm_NEW/chunks/results/region/", isos[i], "_P1_alpha_r_", b, ".tif")))
    alpha10_r[[i]] <- readAll(raster(paste0("calibration/slurm_NEW/chunks/results/region/", isos[i], "_P12_alpha_r_", b, ".tif")))
  } 

#alpha_u_med <- readAll(raster("calibration/slurm_NEW/chunks/results/P12_alpha_u_beta_v4_pos_a_q50.tif"))
#alpha_r_med <- readAll(raster("calibration/slurm_NEW/chunks/results/P12_alpha_r_beta_v4_pos_a_q50.tif"))


# (3) load and preprocess coastal and urban mask data
cst_med <- readAll(raster("boundary_data/coastal/coastal_mask0.tif"))
urb_med <- readAll(raster("boundary_data/urban/ghsl/smod_wgs/urb_mask15.tif"))


# (4) load and preprocess pop data
# -------------------------------------
pop_rur_med <- readAll(stack(paste0("pop_data/pop_", t_pop1, "/pop", t_pop2, "_R2018_rur.tif")))
pop_urb_med <- readAll(stack(paste0("pop_data/pop_", t_pop1, "/pop", t_pop2, "_R2018_urb.tif")))
pop_tot_med <- readAll(stack(paste0("pop_data/pop_", t_pop1, "/pop", t_pop2, "_R2018.tif")))


# (5) load distance raster for reducing the number of grid cells of large countries
# -------------------------------------
# Option 2: produce euclidean distance raster in ArcGIS and load here
buf <- readAll(raster("pop_data/pop_2015/pop15_R2018_eucl.tif")) # distance raster produced in ArcGIS
buf[buf > 30000] <- NA
buf[!is.na(buf)] <- 1


save.image(paste(paste(paste(task, task, sep = "/"), "input", sep = "_"), "RData", sep = "."))

