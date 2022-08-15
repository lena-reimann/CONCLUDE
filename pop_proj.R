##########################################################
## Population projections accounting for sea-level rise ##
##########################################################
# by Lena Reimann
# Aug 15, 2022

# Goal: produce population projections that account for the effects of SLR on spatial changes in population distributions 
# comments: # 1) basic model setup: Reimann et al. (2021) Accounting for internal migration in spatial population projections-a gravity-based modeling approach using the Shared Socioeconomic Pathways. In: Environ Res Lett 16, 074025. DOI: https://doi.org/10.1088/1748-9326/ac0b66
            # 2) SLR-induced migration and adaptation policies model setup: Reimann et al. (2022) xxxxx
            # 3) possible scenario combinations: sp1, ssp3, spp5; nA, wA

rm(list=ls())

lib = "./r_libs"

library("sp", lib.loc = lib)  
library("raster", lib.loc = lib)  
library("geosphere", lib.loc = lib) 
library("scales", lib.loc = lib) 
library("Matrix", lib.loc = lib)
library("rgdal", lib.loc = lib)

#for parallel processing
library("parallel", lib.loc = lib) 
library("foreach", lib.loc = lib) 
library("doParallel", lib.loc = lib) 
library("igraph", lib.loc = lib) 
library("doSNOW", lib.loc = lib) 

#-------------------------------------------------#
#### Step 1 - Define basic information of task ####
#-------------------------------------------------#
task          <- "pop_proj"
id            <- c(8, 12, 70, 818, 376, 422, 434, 499, 504, 275, 760, 788, 792) #ISO3 country codes numeric
iso           <- c("ALB", "DZA", "BIH", "EGY", "ISR", "LBN", "LBY", "MNE", "MAR", "PSE", "SYR", "TUN", "TUR") #ISO3 country codes letters
reg           <- "SE" #geographical region: SE, N            
ssp           <- "ssp3" #SSP to be run: ssp1, ssp3, ssp5             
scen          <- "nA" #adaptation scenario: nA, wA
gw            <- ifelse(reg == "SE", 10, 20) #gravity window size for distance-decay function

cores         <- 5 #define the number of cores to be used for parallel processing (i.e. countries run in parallel)


#--------------------------#
#### Step 2 - Load data ####
#--------------------------#

## define paths for loading the data ##
path <- "path_to_data" # path to input data
path_m <- "path_to_distance_matrix" # path to distance matrix: sparse distance matrix per country of the respective gravity window of the distance-decay effect (see Reimann et al. 2021); produced with the 'Matrix' package


# define time steps
t <- seq(2010, 2100, by=10)  # projections years including 2010 (baseline)
t_proj <- seq(2020,2100, by= 10) # maybe add fove year step later in case the process is fast (but: SLR variables also only produced in 10-year steps)


## load the data  ##
# load country ids based on Gridded Population of the World (GPW)
gpw_ids <- read.csv("cntr_ids.csv")
isos <- gpw_ids$iso
ids <- gpw_ids$gpw_id


## (1) load data used in all scenario combinations (see Reimann et al. 2021)

# spatial mask for the entire Mediterranean region
mask00_med <- raster(paste(path, "sp_mask.tif", sep = "/"))

# Ai as calibrated and preprocessed for urban versus rural locations
alpha_u_med <- raster(paste(path, "alpha_u.tif", sep = "/"))
alpha_r_med <- raster(paste(path, "alpha_r.tif", sep = "/"))

# coastal mask and urban mask
cst_med <- raster(paste(path, "cst.tif", sep = "/")) #defining coastal areas (= Low Elevation Coastal Zone (LECZ) + 20 km coastline buffer)
urb_med <- raster(paste(path, "urb.tif", sep = "/")) #defining urban areas (based on SMOD and adjusted to the country-specific definition of urban population)

# population base year data
pop_rur_med <- raster(paste(path, "pop_rur.tif", sep = "/")) #rural
pop_urb_med <- raster(paste(path, "pop_urb.tif", sep = "/")) #urban
pop_tot_med <- raster(paste(path, "pop_tot.tif", sep = "/")) #total

# distance raster produced with euclidean distance tool for reducing the number of grid cells of large countries in uninhabitable locations
buf <- raster(paste(path, "distance.tif", sep = "/")) 
buf[buf > 30000] <- NA
buf[!is.na(buf)] <- 1

#  load ssp projections file with all ssps, time steps, and countries
ssp_proj <- read.csv(paste(path, "ssp_proj.csv", sep = "/"), stringsAsFactors = F) 


## (2) load data specific to the ssp, rcp, and adaptation scenario to be run (see Reimann et al. 2022)

if (ssp == "ssp1") {
  # define t for spatial mask (loaded further down): leave out the year 2010 (because we start projection for 2020 and the 2010 mask already accounts for SLR)
  t_sp_mask <- c(2000, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100) 
  
  # inundation rasters per ssp-rcp combination and adaptation scenario
  if (scen == "nA") { #submergence per t based on rcp2.6
    name <- paste0(paste(paste(path, "inun_ssp1-rcp26_nA", sep = "/"), t, sep = "_"), ".tif")
    inun <- stack(name)
  } else { #inundation raster of retreat return period per t based on rcp2.6 in unprotected locations (# we do not need the submergence layer because already part of the retreat layer)
    name <- paste0(paste(paste(path, "inun_ssp1-rcp26_wA", sep = "/"), t, sep = "_"), ".tif")
    inun <- stack(name)
  }
  
  } else if (ssp == "ssp3") {
    # define t for spatial mask (loaded further down): 2000 (no change in habitability)
    t_sp_mask <- c(2000) 
    
    # inundation rasters per ssp-rcp combination and adaptation scenario
    if (scen == "nA") { #submergence per t based on rcp4.5
      name <- paste0(paste(paste(path, "inun_ssp3-rcp45_nA", sep = "/"), t, sep = "_"), ".tif")
      inun <- stack(name)
    } else { #submergence per t based on rcp4.5 in unprotected locations
      name <- paste0(paste(paste(path, "inun_ssp3-rcp45_wA", sep = "/"), t, sep = "_"), ".tif")
      inun <- stack(name)
    }
    
  } else {
    #define t for spatial mask (loaded further down): 2000, 2010 (coastline buffer only)
    t_sp_mask <- c(2000, 2010)
    
    # inundation rasters per ssp-rcp combination and adaptation scenario
    if (scen == "nA") { #submergence per t based on rcp8.5
      name <- paste0(paste(paste(path, "inun_ssp5-rcp85_nA", sep = "/"), t, sep = "_"), ".tif")
      inun <- stack(name)
    } else { #submergence per t based on rcp8.5 in unprotected locations
      name <- paste0(paste(paste(path, "inun_ssp5-rcp85_wA", sep = "/"), t, sep = "_"), ".tif")
      inun <- stack(name)
    }
}

# load spatial mask
mask_cntr <- list()
for (i in seq(ids)) {
  name <- paste0(paste(paste0(paste(paste(path, ssp, "sp_mask", sep = "/"), ssp, "id", sep = "_"), ids[i]), t_sp_mask, sep = "_"), ".tif")  
  mask_cntr[[i]] <- stack(name)
}



#---------------------------------#
#### Step 3 - define functions ####
#---------------------------------#

set0 <- function(x) {x[ is.na(x)] <- 0; return(x)}

## I - pop projections ##
pop_proj_SLR <- function(beta_cst, beta_inl, vec_alpha, pop_dens, cells_dist, cst_mask, pop_base_tot,
                     cntr_li, chg, pop_base, inun) {
  # (0) implement inundated area due to SLR
  # 'pick up' population in inundated locations and add to chg
  pop_inun <- pop_base * inun
  pop_tot_inun <- pop_base_tot * inun
  
  pop_base <- pop_base - pop_inun # 'overwrite' pop with pop minus the inundated pop
  pop_base_tot <- pop_base_tot - pop_tot_inun
 
  chg <- chg + sum(pop_inun) # add to chg
  
  #update cntr_li to decrease human habitation potential in inundated locations
  cntr_li <- cntr_li * (1-inun)
  
  
  # (1) calculate potential for the entire country based on cst and inl pot
  # calculate pje cst
  exp_cells_dist <- cells_dist
  exp_cells_dist@x <- exp(-beta_cst * exp_cells_dist@x)
  pje_cst <- as.vector(exp_cells_dist %*% pop_base_tot) 
  
  # calculate pje inl
  exp_cells_dist <- cells_dist
  exp_cells_dist@x <- exp(-beta_inl * exp_cells_dist@x)
  pje_inl <- as.vector(exp_cells_dist %*% pop_base_tot) 
  
  # calculate potentials without Li first
  pot <- ifelse(cst_mask > 0.5, 
         pje_cst + (pop_base_tot * vec_alpha), 
         pje_inl + (pop_base_tot * vec_alpha))

  # condition for using Li: do not use Li in case potential is negative  
  pot <- ifelse(pot < 0 & cntr_li == 0, pot, pot * cntr_li) #new version resulting in a lower likelihood of sum(pot) <0
  #pot <- ifelse(pot > 0, cntr_li * pot, pot) #old and wrong version
  
 
  # (2) assign pop changes to each cell
  if (chg >= 0) { # if pop change positive
    # if sum(pot) below 0, set pot<0 to 0 
    ifelse(sum(pot) < 0, pot[pot < 0] <- 0, pot)
    
    # potential of each cell relative to full potential
    rchange <- pot / sum(pot, na.rm = T) 
    # pop for t+1
    pop_new <- pop_base + chg * rchange
    
  } else { # if pop change negative
    # define basic variables needed for the loop
    pot <- abs(pot) #use abs values of pot to avoid that those cells with neg pot become positive; not needed any more?
    pop_itr <- pop_base
    pop_loss <- chg
    
    # this loop is to assign population decrease to each cell while preventing negative population from appearing
    iter <- 100 # define number of iterations
    for (j in seq_len(iter)) {
      
      if (sum(pot, na.rm = T) > 0) {  # if pot above 0, remove people proportional to the potential in each cell
      r_decr              <- pot / sum(pot, na.rm = T)
      pop_itr             <- pop_itr + pop_loss * r_decr # pop in t+1 (may include negative values; see next line)
      pop_loss            <- sum(pop_itr[which(pop_itr < 0)]) # sum of negative values (i.e. people that need to be deducted again because these led to negative cell values)
      
      pop_itr[which(pop_itr < 0)] <- 0 # pop for next iteration
      pot[pop_itr==0]         <- 0 # pot for next iteration
      
      # pop for next iteration
      if (abs(pop_loss) < 0.000001) { # if remaining pop below value, break the loop
        pop_new  <- pop_itr
        break
      }
      
      } else { # if pot = 0, remove people proportional to the number of people living in each cell
        r_decr              <- pop_itr / sum(pop_itr, na.rm = T)
        pop_itr             <- pop_itr + pop_loss * r_decr
        pop_loss            <- sum(pop_itr[which(pop_itr < 0)])
        
        pop_itr[which(pop_itr < 0)] <- 0 # pop for next iteration
        pot[pop_itr==0]         <- 0 # pot for next iteration
        
        if (abs(pop_loss) == 0) { # if no more people left to be distributed
          pop_new  <- pop_itr
          break
        }
    } # end of else statement for pot = 0
    } # end of j loop
  } # end of else statement for negative pop growth
  
  
  # (3) a) apply population density threshold and b) set low pop values to 0
  # calculate population to be newly distributed
  pop_add <- sum(pop_new[which(pop_new >= pop_dens)] - pop_dens)
  pop_add <- pop_add + sum(pop_new[which(pop_new < 1)])  # should always be positive
  
  #set potential to 0
  pot[which(pop_new >= pop_dens)] <- 0
  pot[which(pop_new < 1)] <- 0  
  
  #calculate rchange
  ifelse(sum(pot) > 0, rchange <- pot / sum(pot), rchange <- 0) # normalize of sum >0, otherwise 0

   # a) set pop to threshold
  pop_new[which(pop_new >= pop_dens)] <- pop_dens
  # b) set low pop values to 0
  pop_new[which(pop_new < 1)] <- 0
  
  # add pop according to pot
  pop_new <- pop_new + pop_add * rchange
 
  return(pop_new)
  
  } # end of function


## II - urban reclassification ##
urb_rcl <- function(urb_mask, tot_proj_tif, urbnp) {
  require(raster, igraph)
  
  # define pop to be reclassified
  urb_proj <- tot_proj_tif
  
  # mask urb mask data to cntr extent
  urb_mask <- crop(urb_mask, urb_proj)
  urb_mask <- mask(urb_mask, urb_proj)
  
  # mask projections with urban mask
  urb_pr <- mask(urb_proj, urb_mask)
  
  # calculate total and relative difference
  urb_tot <- sum(values(urb_pr), na.rm = T) 
  urb_rel <- sum(values(urb_pr), na.rm = T) / urbnp 
  
  if (urb_rel > 1) { # if urban pop overestimated
    
    ## deduct cells
    # find inner boundary cells
    bou8 <- boundaries(urb_pr, type = 'inner', classes = F, directions = 8)
    
    # remove those bou8 that are one clump, i.e. where the entire settlement would 'disappear'
    #make intermediate bou raster
    bou <- bou8
    bou[!is.na(bou)] <- 1
    
    #produce clumps of bou
    cl_bou <- clump(bou, directions = 8)
    
    #calculate zonals
    zon <- zonal(bou8, cl_bou, fun = 'min')
    
    #extract clump ids where min = 1, i.e. entire clump would be removed
    cl_id <- zon[zon[,"min"]==1, "zone"]
    
    if (length(cl_id) != nrow(zon)) {
      #set clumps not to use to NA
      cl_bou[cl_bou %in% cl_id] <- NA
      
      #set bou8 cells to NA based on cl_bou
      bou8[is.na(cl_bou)] <- NA
    }
    
    # select those cells where bou8 == 1
    cells <- which(values(bou8==1))  # cell positions of cells surrounding the urban pop cells
    v_bou <- urb_pr[bou8==1] # pop value of the boundary cells
    
    # write all relevant data into one dataframe 
    df_r <- data.frame(cells, v_bou)
    
    # sort in ascending order based on lowest pop values
    v_bou_sum <- sum(v_bou)
    df_sort <- data.frame(cell = df_r[order(df_r[,"v_bou"], decreasing = F),"cells"], pop = sort(df_r[,"v_bou"], decreasing = F))
    
    # calculate number of people overestimated
    dif <- urb_tot - urbnp
    
    if (dif > 1) {
      # loop through the pop data and substract the pop values until dif <0
      dif_m <- dif  # difference at start
      dif_j <- vector()
      cells_r <- vector()
      pop_r <- vector()
      
      # condition: if one round of boundaries not enough --> run second round of bou
      if (dif >= sum(df_sort[,"pop"])) {
        dif_j <- dif - sum(df_sort[,"pop"])
        dif_m <- dif_j
        cell_sub <- data.frame(cell = df_sort[,"cell"], pop = df_sort[,"pop"], dif_m = rep(dif_m,nrow(df_sort)), dif_j = rep(dif_j,nrow(df_sort)))
        
        if (is.data.frame(cell_sub)) {
          cell_sub <- cell_sub$cell
          urb_pr[cell_sub] <- NA
        }
        
        #set 0s to na
        #urb_pr[urb_pr==0] <- NA
        
        # loop through second round of cells
        bou8 <- boundaries(urb_pr, type = 'inner', classes = F, directions = 8)
        
        # remove those bou8 that are one clump, i.e. where the entire settlement would 'disappear'
        #make intermediate bou raster
        bou <- bou8
        bou[!is.na(bou)] <- 1
        
        #produce clumps of bou
        cl_bou <- clump(bou, directions = 8)
        
        #calculate zonals
        zon <- zonal(bou8, cl_bou, fun = 'min')
        
        #extract clump ids where min = 1, i.e. entire clump would be removed
        cl_id <- zon[zon[,"min"]==1, "zone"]
        
        #set clumps not to use to NA
        cl_bou[cl_bou %in% cl_id] <- NA
        
        #set bou8 cells to NA based on cl_bou
        bou8[is.na(cl_bou)] <- NA
        
        # select those cells where bou8 == 1
        cells <- which(values(bou8==1))  # cell positions of cells surrounding the urban pop cells
        v_bou <- urb_pr[bou8==1] # pop value of the boundary cells
        
        # write all relevant data into one dataframe 
        df_r <- data.frame(cells, v_bou)
        
        # sort in ascending order based on lowest pop values
        v_bou_sum <- sum(v_bou)
        df_sort <- data.frame(cell = df_r[order(df_r[,"v_bou"], decreasing = F),"cells"], pop = sort(df_r[,"v_bou"], decreasing = F))
        
        dif_j <- vector()
        cells_r <- vector()
        pop_r <- vector()
        for (m in 1:nrow(df_sort)) {  
          cell <- df_sort[m,"cell"]
          pop <- df_sort[m,"pop"]
          dif_m <- dif_m - pop
          dif_j <- c(dif_j,dif_m)
          cells_r <- c(cells_r,cell)
          pop_r <- c(pop_r, pop)
          if (dif_m <= 0) { 
            if (m == 1) {
              cell_sub <- data.frame(cells = cells_r, pop = pop_r, dif_m = dif_m, dif_j = dif_j)
            } else if (abs(dif_m) < dif_j[m-1]) {
              # write result into cell_add: col1 = cells to be excluded, col2 = pop, col3 = decreasing difference during loop
              cell_sub <- data.frame(cell = cells_r, pop = pop_r, dif_m = rep(dif_m,length(cells_r)), dif_j = dif_j)
            } else {
              cell_sub <- data.frame(cell = cells_r[1:length(cells_r)-1], pop = pop_r[1:length(pop_r)-1], dif_m = rep(dif_m,(length(cells_r)-1)), 
                                     dif_j = dif_j[1:length(dif_j)-1])
            }
            break
          }
        }
        
        if (is.data.frame(cell_sub)) {
          cell_sub <- cell_sub$cell
          urb_pr[cell_sub] <- NA
        }
        
      } else { # if one round of bou is enough
        # loop through the pop data and substract the pop values until dif <0
        for (m in 1:nrow(df_sort)) {  
          cell <- df_sort[m,"cell"]
          pop <- df_sort[m,"pop"]
          dif_m <- dif_m - pop
          dif_j <- c(dif_j,dif_m)
          cells_r <- c(cells_r,cell)
          pop_r <- c(pop_r, pop)
          if (dif_m <= 0) { 
            if (m == 1) {
              cell_sub <- data.frame(cells = cells_r, pop = pop_r, dif_m = dif_m, dif_j = dif_j)
            } else if (abs(dif_m) < dif_j[m-1]) {
              # write result into cell_add: col1 = cells to be excluded, col2 = pop, col3 = decreasing difference during loop
              cell_sub <- data.frame(cell = cells_r, pop = pop_r, dif_m = rep(dif_m,length(cells_r)), dif_j = dif_j)
            } else {
              cell_sub <- data.frame(cell = cells_r[1:length(cells_r)-1], pop = pop_r[1:length(pop_r)-1], dif_m = rep(dif_m,(length(cells_r)-1)), 
                                     dif_j = dif_j[1:length(dif_j)-1])
            }
            break
          }
        } # end of m loop
        
        if (is.data.frame(cell_sub)) {
          cell_sub <- cell_sub$cell
          urb_pr[cell_sub] <- NA
        }
      }
      
    } # end of if dif > 1 condition
    
  } else if (urb_rel < 1) { # if urban pop underestimated 
    
    # calculate number of people underestimated
    dif <- urbnp - urb_tot
    
    if (dif > 1) {
      
      ## add cells: run loop x times (max 10)
      itr <- 10
      
      for (n in 1:itr) {
        # find outer boundary cells
        bou8 <- boundaries(urb_pr, type = 'outer', classes = F, directions = 8)
        
        # add all pop data where bou == 1 to the raster
        urb_add <- urb_pr
        urb_add[bou8==1] <- urb_proj[bou8 == 1]
        
        # select those cells where bou8 == 1
        cells <- which(values(bou8==1))  # cell positions of cells surrounding the urban pop cells
        v_bou <- urb_add[bou8==1] # pop value of the boundary cells
        v_bou[is.na(v_bou)] <- 0
        
        # write all relevant data into one dataframe 
        df_r <- data.frame(cells, v_bou)
        
        # sort in descending order based on highest pop values
        v_bou_sum <- sum(v_bou)
        df_sort <- data.frame(cell = df_r[order(df_r[,"v_bou"], decreasing = T),"cells"], pop = sort(df_r[,"v_bou"], decreasing = T))
        
        # condition: if one round of boundaries not enough 
        if (dif >= sum(df_sort[,"pop"])) {
          dif_j <- dif - sum(df_sort[,"pop"])
          dif_m <- dif_j
          cell_add <- data.frame(cell = df_sort[,"cell"], pop = df_sort[,"pop"], dif_m = rep(dif_m,nrow(df_sort)), dif_j = rep(dif_j,nrow(df_sort)))
          
          if (is.data.frame(cell_add)) {
            cell_add <- cell_add$cell
            urb_pr[cell_add] <- urb_proj[cell_add]
          }
          
          # calculate number of people still missing -> update dif
          dif <- dif_m      
        } else {
          # loop through the pop data and add the pop values until dif <0
          dif_m <- dif  # difference at start
          dif_j <- vector()
          cells_r <- vector()
          pop_r <- vector()
          
          for (m in 1:nrow(df_sort)) {  
            cell <- df_sort[m,"cell"]
            pop <- df_sort[m,"pop"]
            dif_m <- dif_m - pop
            dif_j <- c(dif_j,dif_m)
            cells_r <- c(cells_r,cell)
            pop_r <- c(pop_r, pop)
            if (dif_m <= 0) { 
              if (m == 1) {
                cell_add <- data.frame(cells = cells_r, pop = pop_r, dif_m = dif_m, dif_j = dif_j)
              } else if (abs(dif_m) < dif_j[m-1]) {
                # write result into cell_add: col1 = cells to be excluded, col2 = pop, col3 = decreasing difference during loop
                cell_add <- data.frame(cell = cells_r, pop = pop_r, dif_m = rep(dif_m,length(cells_r)), dif_j = dif_j)
              } else {
                cell_add <- data.frame(cell = cells_r[1:length(cells_r)-1], pop = pop_r[1:length(pop_r)-1], dif_m = rep(dif_m,(length(cells_r)-1)), 
                                       dif_j = dif_j[1:length(dif_j)-1])
              }
              
              if (is.data.frame(cell_add)) {
                cell_add <- cell_add$cell
                urb_pr[cell_add] <- urb_proj[cell_add]
              }
              
              break
            }
            
          } # end of m loop
          
          break 
        }
      } # end of n loop
      
    } # end of if dif > 1 condition
    
  } else { # if not understimated or overestimated
    urb_pr <- urb_pr
  }
  
  # calculate total urban pop for cross-checking
  urb_tot <- sum(zonal(urb_pr, mask00, fun = 'sum')[,2])  
  urb_rel <- sum(zonal(urb_pr, mask00, fun = 'sum')[,2]) /  urbnp
  
  return(urb_pr)
  
} # end of function



#--------------------------------#
#### Step 4 - preprocess data ####
#--------------------------------#

## betas for each ssp ##
# betas established during calibration process
if (reg == "SE") {
  beta_cr <- 0.0840616795915027 # averaged for P12 and adjusted based on urban sprawl analysis
  beta_cu <- 0.524418163450559
  beta_ir <- 0.0929102774432398
  beta_iu <- 0.709506927021345
} else {
  beta_cr <- 1.8 # averaged for P12 and adjusted based on urban sprawl analysis
  beta_cu <- 0.43641450168164
  beta_ir <- 2.2
  beta_iu <- 0.482352870279708
}

# calculate new betas for the scenario 
if (reg == "SE") { # for southern and eastern Med
  if (ssp == "ssp1") {
    beta_cr <- beta_cr * 3 # established beta * correction value for the SSP --> change!?
    beta_cu <- beta_cu * 3
    beta_ir <- beta_ir * 3
    beta_iu <- beta_iu * 3
  } else if (ssp == "ssp2") {
    beta_cr <- beta_cr * 1.2
    beta_cu <- beta_cu * 1.2
    beta_ir <- beta_ir * 1.2
    beta_iu <- beta_iu * 1.2
  } else if (ssp == "ssp3") {
    beta_cr <- beta_cr * 0.98
    beta_cu <- beta_cu * 0.98
    beta_ir <- beta_ir * 0.98
    beta_iu <- beta_iu * 0.98
  } else if (ssp == "ssp4") {
    beta_cr <- beta_cr * 0.9
    beta_cu <- beta_cu * 0.9
    beta_ir <- beta_ir * 0.9
    beta_iu <- beta_iu * 0.9
  } else {
    beta_cr <- beta_cr * 0.85
    beta_cu <- beta_cu * 0.85
    beta_ir <- beta_ir * 0.85
    beta_iu <- beta_iu * 0.85
  }
} else {   # for northern Med
  if (ssp == "ssp1") {
    beta_cr <- beta_cr * 4
    beta_cu <- beta_cu * 4
    beta_ir <- beta_ir * 4
    beta_iu <- beta_iu * 4
  } else if (ssp == "ssp2") {
    beta_cr <- beta_cr * 1.5
    beta_cu <- beta_cu * 1.5
    beta_ir <- beta_ir * 1.5
    beta_iu <- beta_iu * 1.5
  } else if (ssp == "ssp3") {
    beta_cr <- beta_cr * 0.98
    beta_cu <- beta_cu * 0.98
    beta_ir <- beta_ir * 0.98
    beta_iu <- beta_iu * 0.98
  } else if (ssp == "ssp4") {
    beta_cr <- beta_cr * 0.95
    beta_cu <- beta_cu * 0.95
    beta_ir <- beta_ir * 0.95
    beta_iu <- beta_iu * 0.95
  } else {
    beta_cr <- beta_cr * 0.9
    beta_cu <- beta_cu * 0.9
    beta_ir <- beta_ir * 0.9
    beta_iu <- beta_iu * 0.9
  }
}


## pop density thresholds (#established by analyzing the currently (i.e. 2015) highest pop densities in the two regions) ##
if (reg == "SE") {
  pop_dens_u <- 80291
  pop_dens_r <- 46193
} else {
  pop_dens_u <- 67840
  pop_dens_r <- 10011
}

if (ssp == "ssp1") {
  pop_dens_u  <- pop_dens_u * 1.1 # currently observed density * correction value for the SSP
  pop_dens_r  <- pop_dens_r * 1.1
} else if (ssp == "ssp2") {
  pop_dens_u  <- pop_dens_u * 1.05 # higher pop density than currently observed
  pop_dens_r  <- pop_dens_r * 1.05  
} else if (ssp == "ssp4") {
  pop_dens_u  <- pop_dens_u * 0.95 # lower pop density than currently observed
  pop_dens_r  <- pop_dens_r * 0.95
} else if (ssp == "ssp5") {
  pop_dens_u  <- pop_dens_u * 0.9 # lower pop density than currently observed
  pop_dens_r  <- pop_dens_r * 0.9
} else {
  pop_dens_u  <- pop_dens_u # no change
  pop_dens_r  <- pop_dens_r
}



#--------------------------------------------------------------#
#### Step 5 - Loop through each cntr to produce projections ####
#--------------------------------------------------------------#
start <- Sys.time()
start

cl <- makeCluster(cores)
registerDoParallel(cl, cores = cores)

results <- 
foreach(k = 1:length(id),
        .packages = c("raster", "igraph", "Matrix")) %dopar% {

  #----------------------------------------------#
  # I - preprocess data for countries to project #
  #----------------------------------------------#
  
  ## load distance matrix
  dist_m <- readMM(paste(paste(path_m, iso[k], sep = "/"), "mtx.gz", sep = "."))
  
  ## mask ##
  # extract mask data needed
  mask00 <- mask_cntr[[which(id[k]==ids)]][[1]]
  
  # different approaches per scenario combination
  if (scen == "nA" || scen == "wA" && ssp == "ssp3") {
    ## for very large cntr with limited area available for settlement (i.e. dza, lby, egy) ## 
    #reduce number of cells to be processed based on a buffer layer
    if (id[k] == 12 || id[k] == 434 || id[k] == 818) { # if dza, lby, or egy to be projected
      buf30 <- resample(buf, mask00, method = "ngb")
      mask_t <- mask00 
      
      # if dza reduce the buffer to 20 km to reduce processing time
      if (id[k] == 12) {
        buf20 <- buf30
        for (z in 1:10) {
          bou <- boundaries(buf20, type = 'inner', directions = 8)
          buf20[bou==1] <- NA 
        }
        
        mask_t[is.na(buf20)] <- NA
        mask <- rasterToPoints(mask_t) 
        
      } else { # use 30km buffer for lby and egy
        mask_t[is.na(buf30)] <- NA
        mask <- rasterToPoints(mask_t)        
      }
      
    } else { # for all other countries
      mask <- rasterToPoints(mask00)
    }
    
  } else if (scen == "wA" && ssp == "ssp1") {
    # mask per time step to include setback zones -> converted to matrices
    mask <- list()
    for (i in seq(t)) {
      if (id[k] == 12 || id[k] == 434 || id[k] == 818) {
        buf30 <- resample(buf, mask00, method = "ngb")
        mask_t <- mask_cntr[[which(id[k]==ids)]][[i]]
        
        # if dza reduce the buffer to 20km to reduce processing time
        if (id[k] == 12) {
          buf20 <- buf30
          for (z in 1:10) {
            bou <- boundaries(buf20, type = 'inner', directions = 8)
            buf20[bou==1] <- NA 
          }
          
          mask_t[is.na(buf20)] <- NA
          mask[[i]] <- rasterToPoints(mask_t) 
          
        } else { # use 30km buffer for lby and egy
          mask_t[is.na(buf30)] <- NA
          mask[[i]] <- rasterToPoints(mask_t)        
        }
        
      } else { # for all other countries
        mask[[i]] <- rasterToPoints(mask_cntr[[which(id[k]==ids)]][[i]])
      }
    } # end of i loop

  } else if (scen == "wA" && ssp == "ssp5") {
    # mask of 2010 for the projections (includes setback zones as of ICZM protocol)
    mask <- mask_cntr[[which(id[k]==ids)]][[2]]
    
    # mask -> converted to matrix
    if (id[k] == 12 || id[k] == 434 || id[k] == 818) { 
      buf30 <- resample(buf, mask, method = "ngb")
      mask_t <- mask 
      
      # if dza reduce the buffer to 20 km to reduce processing time
      if (id[k] == 12) {
        buf20 <- buf30
        for (z in 1:10) {
          bou <- boundaries(buf20, type = 'inner', directions = 8)
          buf20[bou==1] <- NA 
        }
        
        mask_t[is.na(buf20)] <- NA
        mask <- rasterToPoints(mask_t) 
        
      } else { # use 30km buffer for lby and egy
        mask_t[is.na(buf30)] <- NA
        mask <- rasterToPoints(mask_t)        
      }
      
    } else { # for all other countries
      mask <- rasterToPoints(mask)
    }
  }
  
  cntr_LL <- mask[,c(1,2)] # get lon/lat based on the spatial mask
 
  
  ## alphas (here: already scaled to q50) ##
  # crop scaled alphas to respective country
  alpha_u <- crop(alpha_u_med, mask00) #original version
  alpha_r <- crop(alpha_r_med, mask00)   
  
  alpha_u <- mask(alpha_u, mask00)
  alpha_r <- mask(alpha_r, mask00)
  
  if (id[k] == 12 || id[k] == 434 || id[k] == 818) { # if dza, lby, or egy to be projected
    if (id[k] == 12) {
      alpha_u[is.na(buf20)] <- NA
      alpha_r[is.na(buf20)] <- NA
    } else {
    alpha_u[is.na(buf30)] <- NA
    alpha_r[is.na(buf30)] <- NA
    }
    alpha_u <- rasterToPoints(alpha_u)[,3]
    alpha_r <- rasterToPoints(alpha_r)[,3]
  } else {
  alpha_u <- rasterToPoints(alpha_u)[,3]
  alpha_r <- rasterToPoints(alpha_r)[,3]
  }

  
  ## coast ##
  # crop coastal mask data to respective country
  cst <- crop(cst_med, mask00)
  cst <- mask(cst, mask00)
  
  if (id[k] == 12 || id[k] == 434 || id[k] == 818) { # if dza or lby to be projected
    if (id[k] == 12) {
      cst[is.na(buf20)] <- NA
    } else {
      cst[is.na(buf30)] <- NA
    }
    cst <- rasterToPoints(cst)[,3]
    
    } else {
    cst <- rasterToPoints(cst)[,3]
    }
  
  
  ## urban mask data to respective country
  urb <- crop(urb_med, mask00)

  # mask
  urb <- mask(urb, mask00)

  if (id[k] == 12 || id[k] == 434 || id[k] == 818) { # if dza or lby to be projected
    if (id[k] == 12) {
      urb[is.na(buf20)] <- NA
    } else {
      urb[is.na(buf30)] <- NA
    }
  }
  
  
  ## projections ##
  # extract projections of country used
  proj_cntr <- ssp_proj[which(ssp_proj$iso==iso[k]),1:ncol(ssp_proj)] # change accordingly if other start date
  
  # extract scenario used
  proj_cntr <- proj_cntr[which(proj_cntr$method==ssp),2:ncol(proj_cntr)]
  
  # extract t used no because we need the pop change between 2010 and 2020 as well (but what if we use pop 2010 as starting point?)
  #proj_cntr <- proj_cntr[which(proj_cntr$year %in% t_proj),]
  
  # calculate the pop change per time step for R,U,TOT
  proj_cntr$totchng <- c(proj_cntr[1:nrow(proj_cntr), 3]) - c(NA, proj_cntr[1:nrow(proj_cntr)-1,3])
  proj_cntr$urbchng <- c(proj_cntr[1:nrow(proj_cntr), 5]) - c(NA, proj_cntr[1:nrow(proj_cntr)-1,5])
  proj_cntr$rurchng <- c(proj_cntr[1:nrow(proj_cntr), 6]) - c(NA, proj_cntr[1:nrow(proj_cntr)-1,6])
  
  
  ## pop ##
  # crop pop raster to mask00
  pop_tot <- crop(pop_tot_med, mask00)
  pop_urb <- crop(pop_urb_med, mask00)
  pop_rur <- crop(pop_rur_med, mask00)
  
  # align rasters so that number of cells align
  pop_tot <- resample(pop_tot, mask00, method= "ngb")
  pop_urb <- resample(pop_urb, mask00, method= "ngb")
  pop_rur <- resample(pop_rur, mask00, method= "ngb")
  
  pop_tot[is.na(pop_tot)] <- 0
  pop_urb[is.na(pop_urb)] <- 0
  pop_rur[is.na(pop_rur)] <- 0
  
  pop_tot <- mask(pop_tot, mask00)
  pop_urb <- mask(pop_urb, mask00)
  pop_rur <- mask(pop_rur, mask00)
  
  # convert to points
  if (id[k] == 12 || id[k] == 434 || id[k] == 818) { # if dza or lby to be projected
   if (id[k] == 12) {
     pop_tot[is.na(buf20)] <- NA
     pop_urb[is.na(buf20)] <- NA
     pop_rur[is.na(buf20)] <- NA
   } else {
    pop_tot[is.na(buf30)] <- NA
    pop_urb[is.na(buf30)] <- NA
    pop_rur[is.na(buf30)] <- NA
   }
    
    pop_tot <- as.data.frame(rasterToPoints(pop_tot))
    pop_urb <- as.data.frame(rasterToPoints(pop_urb))
    pop_rur <- as.data.frame(rasterToPoints(pop_rur))
  } else {
    pop_tot <- as.data.frame(rasterToPoints(pop_tot))
    pop_urb <- as.data.frame(rasterToPoints(pop_urb))
    pop_rur <- as.data.frame(rasterToPoints(pop_rur))
  }
  
  
  ## scale pop data to harmonize with base year data of ssps ##
  # total pop
  sc_val <- ifelse(sum(pop_tot[,3]) == 0, 0, proj_cntr[which(proj_cntr$year == t[1]),"totlp"] / sum(pop_tot[,3]))
  pop_tot[,3] <- pop_tot[,3] * sc_val
  
  # urban pop
  sc_val <- ifelse(sum(pop_urb[,3]) == 0, 0, proj_cntr[which(proj_cntr$year == t[1]),"urbnp"] / sum(pop_urb[,3]))
  pop_urb[,3] <- pop_urb[,3] * sc_val
  
  # rural pop
  sc_val <- ifelse(sum(pop_rur[,3]) == 0, 0, proj_cntr[which(proj_cntr$year == t[1]),"rurlp"] / sum(pop_rur[,3]))
  pop_rur[,3] <- pop_rur[,3] * sc_val
  
  
  ## inundation due to SLR ##
  inun_poi <- list()
  # preprocess inun data
  for (i in seq(t)) {
    rast <- crop(inun[[i]], mask00)
    rast[is.na(rast)] <- 0
    rast <- extend(rast, mask00, value = 0) # for countries that extend further than the inun layer
    rast <- mask(rast, mask00)
    
    # convert to points
    if (id[k] == 12 || id[k] == 434 || id[k] == 818) { # if dza or lby to be projected
      if (id[k] == 12) {
        rast[is.na(buf20)] <- NA
      } else {
        rast[is.na(buf30)] <- NA
      }
      #convert to points and extract
      inun_poi[[i]] <- rasterToPoints(rast)[,3] 
    } else {
      #convert to points and extract
      inun_poi[[i]] <- rasterToPoints(rast)[,3] 
    }
    #convert to points and extract
    inun_poi[[i]] <- rasterToPoints(rast)[,3] 
  } # end of i loop
  
  
  # define additional variables needed for the projections
  BR      <- dim(cntr_LL)[1]
  M       <- dim(proj_cntr)[1] # 10 (could also be 11; depends on the csv file created, i.e. the projection steps)
  
  ## create arrays of each time step for 2010-2100 ##
  tot_pop  <- array(0, dim = c(BR, 1, M)) # total population, 230387 x 1 x 10
  urb_pop  <- array(0, dim = c(BR, 1, M)) # urban population, 230387 x 1 x 10
  rur_pop  <- array(0, dim = c(BR, 1, M)) # rural population, 230387 x 1 x 10

  # write the data of the base year into the projections array
  tot_pop[,,1] <- pop_tot[,3]
  urb_pop[,,1] <- pop_urb[,3]
  rur_pop[,,1] <- pop_rur[,3]
  
  #----------------------------------------#
  # II - produce projections for each cntr #
  #----------------------------------------#
  
  
  ### a) produce pop projection ###
  for (m in 2:M) {   
  # differentiate ssp1 wA from the others because of the changing spatial mask
    if (ssp == "ssp1" && scen = "wA") {
      urb_pop[,,m] <- pop_proj_SLR(beta_cu, 
                                   beta_iu, 
                                   alpha_u, 
                                   pop_dens_u, 
                                   dist_m, 
                                   cst, 
                                   tot_pop[,,m-1],
                                   mask[[m]][,3], 
                                   proj_cntr$urbchng[m], 
                                   urb_pop[,,m-1],
                                   inun_poi[[m]]) 
      
      rur_pop[,,m] <- pop_proj_SLR(beta_cr, 
                                   beta_ir, 
                                   alpha_r, 
                                   pop_dens_r, 
                                   dist_m, 
                                   cst, 
                                   tot_pop[,,m-1],
                                   mask[[m]][,3], 
                                   proj_cntr$rurchng[m], 
                                   rur_pop[,,m-1],
                                   inun_poi[[m]]) 
    } else {
       urb_pop[,,m] <- pop_proj_SLR(beta_cu, 
                                    beta_iu, 
                                    alpha_u, 
                                    pop_dens_u, 
                                    dist_m, 
                                    cst, 
                                    tot_pop[,,m-1],
                                    mask[,3], 
                                    proj_cntr$urbchng[m], 
                                    urb_pop[,,m-1],
                                    inun_poi[[m]]) 
    
      rur_pop[,,m] <- pop_proj_SLR(beta_cr, 
                                   beta_ir, 
                                   alpha_r, 
                                   pop_dens_r, 
                                   dist_m, 
                                   cst, 
                                   tot_pop[,,m-1],
                                   mask[,3], 
                                   proj_cntr$rurchng[m], 
                                   rur_pop[,,m-1],
                                   inun_poi[[m]])      
    }

    # add up the rural and urban population to get total population
    tot_pop[,,m] <- urb_pop[,,m] + rur_pop[,,m]
    
    
    ### b) reclassify urb ###
    # make raster 
    tot_tif <- rasterize(cntr_LL, mask00, tot_pop[,,m])
    
    # run urb_rcl function
    urb_pr <- urb_rcl(urb, # assuming urb is a raster
                      tot_tif, 
                      proj_cntr$urbnp[m])
    
    # make urban mask based on the urb pop data
    urb <- urb_pr
    urb[!is.na(urb)] <- 1
    
    
    ## make points and write into arrays again
    urb_pr <- calc(urb_pr, set0)
    urb_pr <- mask(urb_pr, mask00)
    
    if (id[k] == 12 || id[k] == 434 || id[k] == 818) { # if dza or lby to be projected
      if (id[k] == 12) {
        urb_pr[is.na(buf20)] <- NA
      } else {
        urb_pr[is.na(buf30)] <- NA
      }
    }
    
    urb_pop[,,m] <- rasterToPoints(urb_pr, na.rm = T)[,3]
    rur_pop[,,m] <- rasterToPoints(tot_tif)[,3] - urb_pop[,,m]

    
    #test
    #sum(urb_pop[,,m] + rur_pop [,,m]) - sum(tot_pop[,,m])
    
  } # end of m loop
  
  #----------------------#
  # III - write the data #
  #----------------------#
  
  # define model period
  model_period <- seq(proj_cntr[which(proj_cntr[,1] == t_proj[1]),1], proj_cntr[M,1], by = 10)
  
  tot_proj <- pop_tot
  urb_proj <- pop_urb
  rur_proj <- pop_rur
  
  # make df with pop proj based on calibration pop data
  for (i in 1:length(model_period)) {
    #tot
    new_col              <- paste0("tot_pop", model_period[i])
    tot_proj[, new_col]  <- tot_pop[,,i+1]      # create new columns, these new columns are total population from 2010 to 2100.
    
    #urb
    new_col              <- paste0("urb_pop", model_period[i])
    urb_proj[, new_col]  <- urb_pop[,,i+1]      # create new columns, these new columns are total population from 2010 to 2100.
    
    #rur
    new_col              <- paste0("rur_pop", model_period[i])
    rur_proj[, new_col]  <- rur_pop[,,i+1]      # create new columns, these new columns are total population from 2010 to 2100
  }
  
  # write csvs
  write.csv(tot_proj, paste(paste(paste(path, task, ssp, scen, isos[which(iso[k]==isos)], sep = "/"), 
                                  "pop_tot", ssp, scen, sep = "_"), "csv", sep = "."), row.names = F)
  write.csv(urb_proj, paste(paste(paste(path, task, ssp, scen, isos[which(iso[k]==isos)], sep = "/"), 
                                  "pop_urb", ssp, scen,  sep = "_"), "csv", sep = "."), row.names = F)
  write.csv(rur_proj, paste(paste(paste(path, task, ssp, scen, isos[which(iso[k]==isos)], sep = "/"), 
                                  "pop_rur", ssp, scen,  sep = "_"), "csv", sep = "."), row.names = F)
  
  # make rasters again 
  tot_proj_tif <- rasterize(tot_proj[,c(1,2)], mask00, tot_proj[,3:ncol(tot_proj)])
  urb_proj_tif <- rasterize(urb_proj[,c(1,2)], mask00, urb_proj[,3:ncol(urb_proj)])
  rur_proj_tif <- rasterize(rur_proj[,c(1,2)], mask00, rur_proj[,3:ncol(rur_proj)])
  
  # write rasters
  writeRaster(tot_proj_tif, filename=paste(paste(paste(path, task, ssp, scen, isos[which(iso[k]==isos)], sep = "/"), "pop_tot", ssp, scen, t, 
                                                 sep = "_"), "tif", sep = "."), format="GTiff", bylayer = T, overwrite=TRUE)
  writeRaster(urb_proj_tif, filename=paste(paste(paste(path, task, ssp, scen, isos[which(iso[k]==isos)], sep = "/"), "pop_urb", ssp, scen, t, 
                                                 sep = "_"), "tif", sep = "."), format="GTiff", bylayer = T, overwrite=TRUE)
  writeRaster(rur_proj_tif, filename=paste(paste(paste(path, task, ssp, scen, isos[which(iso[k]==isos)], sep = "/"), "pop_rur", ssp, scen, t, 
                                                 sep = "_"), "tif", sep = "."), format="GTiff", bylayer = T, overwrite=TRUE)
  
  stack(tot_proj_tif, urb_proj_tif, rur_proj_tif)
  
} #end of k loop

stopCluster(cl) 

end <- Sys.time()
print(end - start)


# ----------------------------- #
#### Step 6 - mosaic rasters ####
# ----------------------------- #

tot_proj_tif <- list()
urb_proj_tif <- list()
rur_proj_tif <- list()

for (i in 1:length(results)) {
  tot_proj_tif[[i]] <- results[[i]][[ 1:10]]
  urb_proj_tif[[i]] <- results[[i]][[11:20]]
  rur_proj_tif[[i]] <- results[[i]][[21:30]]
}

# make mosaic
tot_proj_tif$fun <- sum
urb_proj_tif$fun <- sum
rur_proj_tif$fun <- sum

tot_proj_tif$na.rm <- TRUE
urb_proj_tif$na.rm <- TRUE
rur_proj_tif$na.rm <- TRUE

# make one raster for the entire chunk
tot_reg <- do.call(mosaic, tot_proj_tif)
urb_reg <- do.call(mosaic, urb_proj_tif)
rur_reg <- do.call(mosaic, rur_proj_tif)

# write rasters as grd
writeRaster(tot_reg, filename=paste(paste(paste(path, task, ssp, scen, reg, sep = "/"), "pop_tot", ssp, scen, t, sep = "_"), "tif", 
                                    sep = "."), format="GTiff", bylayer = T, overwrite=TRUE)
writeRaster(urb_reg, filename=paste(paste(paste(path, task, ssp, scen, reg, sep = "/"), "pop_urb", ssp, scen, t, sep = "_"), "tif", 
                                    sep = "."), format="GTiff", bylayer = T, overwrite=TRUE)
writeRaster(rur_reg, filename=paste(paste(paste(path, task, ssp, scen, reg, sep = "/"), "pop_rur", ssp, scen, t, sep = "_"), "tif", 
                                    sep = "."), format="GTiff", bylayer = T, overwrite=TRUE)


# save entire workspace
save.image(paste(paste(paste(path, task, ssp, scen, task, sep = "/"), reg, ssp, scen, sep = "_"), "RData", sep = "."))

