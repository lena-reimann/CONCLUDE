#######################################
## Population projections input data ##
#######################################
# by Lena Reimann
# Nov 9, 2020

# Goal: preprocess input data for population projections for each ssp individually
# Input data: betas (cu, cr, iu, ir), pop data (whole region differentiation between N and SE?)
# comments: # 1) for calculation on rz cluster: two scripts --> (1) data prep; (2) parallel processing
            # 2) v1: with ghsl split into u and r
            # 3) for entire region
            # 4) use readAll() function for including all raster info in the workspace when saving it 
                 #(so that it can be used on another machine and/or on the RZ cluster)
         

rm(list=ls())

#load packages
#library(rgdal)
library(raster)
#library(tmaptools)
#library(geosphere)

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
task          <- "pop_proj"
ssp           <- "ssp4"     # change accordingly


path          <- "D:/PhD/Fulbright/Demographic_model/model_extension" # change accordingly
setwd(paste(path, sep = "/"))       # may keep this line unchanged (set working directory) 
getwd()

# define time steps
t_pop1 <- c(1990, 2000, 2010, 2015)
t_pop2 <-  c("90", "00", "10", "15")
t_pop3 <-  c("pop90", "pop00", "pop10", "pop15")

t <- seq(2010, 2100, by=10)  # projections years including 2010 (baseline)
t1 <- seq(2000, 2100, by=10) # projections years including 2000 (baseline used for model calibration)
t2 <- c(2000, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100) # for SSP1, leaving out the mask of 2010 because it has been produced based on SLR data already; the mask of 2000 is representative of the 2010 mask without SLR
 
t_proj <- seq(2020,2100, by= 10) # maybe add fove year step later in case the process is fast (but: SLR variables also only produced in 10-year steps)


# define gpw ids
gpw_ids <- read.csv("boundary_data/gpw4_10/gpw_id_med_iso25.csv")
isos <- gpw_ids$iso
ids <- gpw_ids$gpw_id


#--------------------------------------------------------------
# Section II  --  get the data
#--------------------------------------------------------------
# (1) load mask data
# for ssp4: load 2000 mask based on ssp5

#entire med
mask00_med <- readAll(raster("spatial_mask/final/sp_mask_ssp5v5.tif"))

#per cntr
mask_cntr <- list()
t_sp_mask <- c(2000)


#load masks per country for each t based on ssp5 (because similar technological davelopments)
for (i in seq(ids)) {
  name <- paste0(paste(paste0(paste(paste0("spatial_mask/sp_mask_cntr25_ssp5/sp_mask"), "ssp5_id", sep = "_"), ids[i]), 
                       t_sp_mask, sep = "_"), ".tif")  
  mask_cntr[[i]] <- readAll(stack(name))
}


# (2) load and preprocess alphas (for now load ghsl of 2015 only, but also add esl data)
alpha_u_med <- readAll(raster("calibration/slurm_NEW/chunks/results/P12_alpha_u_beta_v4_pos_a_q50.tif"))
alpha_r_med <- readAll(raster("calibration/slurm_NEW/chunks/results/P12_alpha_r_beta_v4_pos_a_q50.tif"))


# (3) load and preprocess coastal and urban mask data
cst_med <- readAll(raster("boundary_data/coastal/coastal_mask0.tif"))
urb_med <- readAll(raster("boundary_data/urban/ghsl/smod_wgs/urb_mask15.tif"))


# (4) load ssp projections file
# ------------------------------------------
ssp_proj <- read.csv("pop_data/ssps/proj_model_input25.csv", stringsAsFactors = F)


# (5) load and preprocess pop data
# -------------------------------------
pop_rur_med <- stack(paste0("pop_data/pop_", t_pop1, "/pop", t_pop2, "_R2018_rur.tif"))
pop_urb_med <- stack(paste0("pop_data/pop_", t_pop1, "/pop", t_pop2, "_R2018_urb.tif"))
pop_tot_med <- stack(paste0("pop_data/pop_", t_pop1, "/pop", t_pop2, "_R2018.tif"))

# select base year data here: 2010
pop_rur_med <- readAll(pop_rur_med[[3]])
pop_urb_med <- readAll(pop_urb_med[[3]])
pop_tot_med <- readAll(pop_tot_med[[3]])


# (5b) load distance raster for reducing the number of grid cells of large countries (i.e. dza, lby)
# -------------------------------------
# Option 2: produce euclidean distance raster in ArcGIS and load here
buf <- readAll(raster("pop_data/pop_2015/pop15_R2018_eucl.tif")) # distance raster produced in ArcGIS
buf[buf > 30000] <- NA
buf[!is.na(buf)] <- 1


# (6) load inundation layers for each time step for locations permanently inundated due to SLR (not here because we don't use ssp4 for migration paper)
#--------------------------------------------------------------
#inun_nA <- list()
#inun_wA <- list()

#for (i in seq(t2)) {  
  ## nA
  #load inun for each t (years 2010 to 2100)
#  name <- paste0(paste(paste0("inundation/results/inun30_sl/inun30_sl85"), t[i], sep = "_"), ".tif")
#  inun_nA[[i]] <- readAll(raster(name))
  
  ## wA
  #load inun for each t (years 2010 to 2100)
#  name <- paste0(paste(paste0("inundation/results/inun30_sl_wA/inun30_sl85_wA"), t[i], sep = "_"), ".tif")
#  inun_wA[[i]] <- readAll(raster(name))
#}

save.image(paste(paste(paste(task, task, sep = "/"), "input", ssp,  sep = "_"), "RData", sep = "."))


