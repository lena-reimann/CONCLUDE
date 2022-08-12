#################################
## model validation for chunks ##
#################################
# by Lena Reimann
# Jan 5, 2021

# Goal: validate model for 2010 and 2015 based on calibrated betas and alphas
# comments: # 1) q50: alphas scaled to the interquartile range (i.e. middle 50 %/ Q3-Q1)

rm(list=ls())

.libPaths("~/r_libs")

library("sp", lib.loc = "~/r_libs")  
library("raster", lib.loc = "~/r_libs")  
library("scales", lib.loc = "~/r_libs") 
library("Matrix", lib.loc = "~/r_libs")
library("rgdal", lib.loc = "~/r_libs")

#for parallel processing
library("parallel", lib.loc = "~/r_libs") 
library("foreach", lib.loc = "~/r_libs") 
library("doParallel", lib.loc = "~/r_libs") 
library("igraph", lib.loc = "~/r_libs") 
library("doSNOW", lib.loc = "~/r_libs") 

#-------------------------------------------#
# Step 1 - Define basic information of task #
#-------------------------------------------#
task          <- "validation"
id            <- c(8, 12, 70, 818, 376, 422, 434, 499, 504, 275, 760, 788, 792)
iso           <- c("ALB", "DZA", "BIH", "EGY", "ISR", "LBN", "LBY", "MNE", "MAR", "PSE", "SYR", "TUN", "TUR")
reg           <- "SE"                
ssp           <- "ssp3"              
scen          <- "nSLR"                
model_period  <- 2010 # change accordingly
p             <- "p1"
a             <- "a_q50"

gw            <- ifelse(reg == "SE", 10, 20)


## load the data based on where the script is run ##
run <- "rz"    # change accordingly 

if (run == "rz") {
  cores <- 7      # we can use up to 12!
  path <- "/work_beegfs/sungg688"
  
  load(paste(paste(paste(path, task, task, sep = "/"), "input", sep = "_"), "RData", sep = "."))
  
  path <- "/work_beegfs/sungg688" # set path again
  path_m <- "/work_beegfs/sungg688/dist_m" # path for distance matrix
} else {
  cores <- 4
  #path          <- "E:/model_extension"  #for laptop
  path <- "D:/PhD/Fulbright/Demographic_model/model_extension" # for PC/hardrive
  setwd(paste(path, sep = "/"))
  
  load(paste(paste(paste(path, task, task, sep = "/"), "input", sep = "_"), "RData", sep = "."))
  
  path <- "D:/PhD/Fulbright/Demographic_model/model_extension"
  path_m <- "Theo/dist_m"
}


#---------------------------------#
#### Step 2 - define functions ####
#---------------------------------#

set0 <- function(x) {x[ is.na(x)] <- 0; return(x)}

## 0 - scale alphas ## (unscaled alphas loaded here)
# alpha = vector
# quantiles = text
scale_alpha <- function(alpha, quantiles) { 
  
  # 1. calculate quantiles
  qu_neg <- quantile(alpha[alpha<0], probs = c(.05, .10, .25, .33, .50, .66, .75, .83, .90, .95))
  qu_pos <- quantile(alpha[alpha>0], probs = c(.05, .10, .25, .33, .50, .66, .75, .83, .90, .95))
  
  if (quantiles == "a_q50") {
    # remove upper and lower 25%
    alpha[alpha > qu_pos["75%"]] <- qu_pos["75%"] 
    alpha[alpha < qu_neg["25%"]] <- qu_neg["25%"] 
  } else {
    alpha[alpha > qu_pos["90%"]] <- qu_pos["90%"] 
    alpha[alpha < qu_neg["10%"]] <- qu_neg["10%"] 
  }
  
  # 2. scale alphas from 0-100
  max <- max(alpha)
  min <- min(alpha)
  # establish "max" for scaling factor
  max <- ifelse(max > abs(min), max, abs(min))
  # calc scaling factor
  factor <- ifelse(max == 0, 0, 100 / max)
  # scale alphas (preserve skewness of the function)
  alpha_sc <- alpha * factor
  
  return(alpha_sc)
}


## I - pop projections ##

pop_proj_nSLR <- function(beta_cst, beta_inl, vec_alpha, pop_dens, cells_dist, cst_mask, pop_base_tot,
                     cntr_li, chg, pop_base) {

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
  pot <- ifelse(pot < 0 & cntr_li == 0, pot, pot * cntr_li)
 
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



## III - combine outputs of foreach loop in two lists ##
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


## IV - calculate error terms ##
# a) vector-based 
# m = modeled 
# o = observed
error <- function(m,o) {
  abs <- m - o # total absolute error?
  abs_sum <- sum(abs(abs))
  abs_mn <- abs_sum / length(o) # absolute error divided by the total number of cells (-> mean absolute error)
  rtae <- abs_sum / sum(o) # relative total absolute error
  wmape <- sum(ifelse(o > 0, abs(abs) / o, abs(abs) / 1) * o) / sum(o) # wmape (= weighted mean absolute percentage error)
  
  rmse <- sqrt(mean((m - o)^2))
  rmse_p <- (sqrt(mean((m - o)^2))) / (mean(o)) #rmse%
  rsq <-  1 - ((sum((o - m)^2)) / (sum((o - mean(o))^2)))
  
  return(c(abs_sum, abs_mn, rtae, wmape, rmse, rmse_p, rsq))
}


# ------------------------------ #
#### Step 3 - preprocess data ####
# ------------------------------ #

# calculate new betas for the scenario 
# make betas of P1 and P12
# see beta cali insights for more info on betas used; adjusted based on urban sprawl analysis
if (reg == "SE") {
  if (model_period == 2010) { #v4
    beta_cr <- 0.004875
    beta_cu <- 0.543511617
    beta_ir <- 0.005125
    beta_iu <- 0.729270887
  } else {
    beta_cr <- 0.0840616795915027 # averaged for P12
    beta_cu <- 0.524418163450559
    beta_ir <- 0.0929102774432398
    beta_iu <- 0.709506927021345
  }
} else {
  if (model_period == 2010) {
    beta_cr <- 1.82 # calibrated for P1
    beta_cu <- 0.344505198
    beta_ir <- 2.18
    beta_iu <- 0.380768903
  } else {
    beta_cr <- 1.8 # averaged for P12 and adjusted based on urban sprawl analysis
    beta_cu <- 0.43641450168164
    beta_ir <- 2.2
    beta_iu <- 0.482352870279708
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


## select alphas
if (model_period == 2010) {
  alpha_u_l <- alpha00_u
  alpha_r_l <- alpha00_r
} else {
  alpha_u_l <- alpha10_u
  alpha_r_l <- alpha10_r
}


## select pop
if (model_period == 2010) {
  pop_tot_med <- pop_tot_med[[2:3]]
  pop_urb_med <- pop_urb_med[[2:3]]
  pop_rur_med <- pop_rur_med[[2:3]]
  
  t <- c(2000, 2010)
} else {
  pop_tot_med <- pop_tot_med[[3:4]]
  pop_urb_med <- pop_urb_med[[3:4]]
  pop_rur_med <- pop_rur_med[[3:4]]
  
  t <- c(2010, 2015)
}


#--------------------------------------------------------------#
#### Step 4 - Loop through each cntr to produce projections ####
#--------------------------------------------------------------#
start <- Sys.time()
start

cl <- makeCluster(cores)
registerDoParallel(cl, cores = cores)

results <- 
foreach(k = 1:length(id),
        .combine = 'comb', .multicombine = T, .init = list(list(), list()), # for combining two lists in output
        .packages = c("raster", "igraph", "Matrix")) %dopar% {

  #----------------------------------------------#
  # I - preprocess data for countries to project #
  #----------------------------------------------#
  
  ## load distance matrix
  dist_m <- readMM(paste(paste(path_m, iso[k], sep = "/"), "mtx.gz", sep = "."))
  
  ## mask ##
  # extract mask data needed
  mask00 <- mask_cntr[[which(id[k]==ids)]][[1]]
  
  ## for very large cntr with limited area available for settlement (i.e. dza, lby) ## 
  #reduce number of cells to be processed based on a buffer layer
  if (id[k] == 12 || id[k] == 434 || id[k] == 818) { # if dza or lby to be projected
    buf30 <- resample(buf, mask00, method = "ngb")
    mask <- mask00
    
    # if dza reduce the buffer to 25 km to reduce processing time
    if (id[k] == 12) {
      buf20 <- buf30
      for (i in 1:10) {
        bou <- boundaries(buf20, type = 'inner', directions = 8)
        buf20[bou==1] <- NA 
      }
    
      mask[is.na(buf20)] <- NA
      mask <- rasterToPoints(mask) 
      
    } else { # use 30km buffer for lby and egy
      mask[is.na(buf30)] <- NA
      mask <- rasterToPoints(mask)        
    }
    
  } else { # for all other countries
    mask <- rasterToPoints(mask00)
  }
 
  cntr_LL <- mask[,c(1,2)] # get lon/lat based on the spatial mask
  cntr_Li <- mask[,3]  # get spatial mask values
  
  
  ## alphas ##
  alpha_u <- alpha_u_l[[which(id[k]==ids)]]
  alpha_r <- alpha_r_l[[which(id[k]==ids)]]
  
  #alpha_u <- crop(alpha_u, mask00) #original version
  #alpha_r <- crop(alpha_r, mask00)   
  
  #alpha_u <- mask(alpha_u, mask00)
  #alpha_r <- mask(alpha_r, mask00)
  
  if (id[k] == 12 || id[k] == 434 || id[k] == 818) { # if dza or lby to be projected
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
  
  # scale alphas
  alpha_u <- scale_alpha(alpha_u, a)
  alpha_r <- scale_alpha(alpha_r, a)
  
  ####---------------- for testing set alphas to 1 -------------------####
  #alpha_u <- rep(1, length(alpha_u))
  #alpha_r <- rep(1, length(alpha_r))
  
  
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
  
  ## calculate pop change to be projected ##
  tot_chg <- sum(pop_tot[,4] - pop_tot[,3])
  urb_chg <- sum(pop_urb[,4] - pop_urb[,3])
  rur_chg <- sum(pop_rur[,4] - pop_rur[,3])
  # urb_chg + rur_chg - tot_chg
  
  # calculate urban pop in t+1
  urbnp <- sum(pop_urb[,4])
  
  # define additional variables needed for the projections
  BR      <- dim(cntr_LL)[1]
  M       <- 2 # 10 (could also be 11; depends on the csv file created, i.e. the projection steps)
  
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
  
  for (m in 2:M) { 
    ### a) produce pop projection ###
    urb_pop[,,m] <- pop_proj_nSLR(beta_cu, 
                                  beta_iu, 
                                  alpha_u, 
                                  pop_dens_u, 
                                  dist_m, 
                                  cst, 
                                  tot_pop[,,m-1],
                                  cntr_Li, 
                                  urb_chg, 
                                  urb_pop[,,m-1]) 
    
    rur_pop[,,m] <- pop_proj_nSLR(beta_cr, 
                                  beta_ir, 
                                  alpha_r, 
                                  pop_dens_r, 
                                  dist_m, 
                                  cst, 
                                  tot_pop[,,m-1],
                                  cntr_Li, 
                                  rur_chg, 
                                  rur_pop[,,m-1]) 
    
    # add up the rural and urban population to get total population
    tot_pop[,,m] <- urb_pop[,,m] + rur_pop[,,m]
    
    
    ### b) reclassify urb ###
    # make raster 
    tot_tif <- rasterize(cntr_LL, mask00, tot_pop[,,m])
    
    # run urb_rcl function
    urb_pr <- urb_rcl(urb,
                      tot_tif, 
                      urbnp)
    
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
    
    #---------------------------------------#
    # IV - calculate error terms (per cntr) #
    #---------------------------------------# 
    
    err_tot <- error(tot_pop[,,m], tot_pop[,,m-1])
    err_urb <- error(urb_pop[,,m], urb_pop[,,m-1])
    err_rur <- error(rur_pop[,,m], rur_pop[,,m-1])
    
    err <- rbind(err_tot, err_urb, err_rur)
    
  } # end of m loop
  
  #----------------------#
  # III - write the data #
  #----------------------#
  
  tot_proj <- pop_tot[,c(1:3)]
  urb_proj <- pop_urb[,c(1:3)]
  rur_proj <- pop_rur[,c(1:3)]
  
  # make df
  #tot
  new_col              <- paste0("tot_pop", model_period)
  tot_proj[, new_col]  <- tot_pop[,,2]
    
  #urb
  new_col              <- paste0("urb_pop", model_period)
  urb_proj[, new_col]  <- urb_pop[,,2] 
    
  #rur
  new_col              <- paste0("rur_pop", model_period)
  rur_proj[, new_col]  <- rur_pop[,,2]
    
  # make rasters again 
  tot_proj_tif <- rasterize(tot_proj[,c(1,2)], mask00, tot_proj[,3:ncol(tot_proj)])
  urb_proj_tif <- rasterize(urb_proj[,c(1,2)], mask00, urb_proj[,3:ncol(urb_proj)])
  rur_proj_tif <- rasterize(rur_proj[,c(1,2)], mask00, rur_proj[,3:ncol(rur_proj)])
  
  if (id[k] == 12 || id[k] == 434 || id[k] == 818) { 
    # set NAs to 0
    tot_proj_tif <- calc(tot_proj_tif, set0)
    urb_proj_tif <- calc(urb_proj_tif, set0)
    rur_proj_tif <- calc(rur_proj_tif, set0)
    
    # mask with original spatial mask
    tot_proj_tif <- mask(tot_proj_tif, mask00)
    urb_proj_tif <- mask(urb_proj_tif, mask00)
    rur_proj_tif <- mask(rur_proj_tif, mask00)
  }
  
  # write rasters
  #writeRaster(tot_proj_tif, filename=paste(paste(paste(path, task, ssp, scen, isos[which(iso[k]==isos)], sep = "/"), "pop_tot", ssp, scen, t, p,
  #                                               sep = "_"), "tif", sep = "."), format="GTiff", bylayer = T, overwrite=TRUE)
  #writeRaster(urb_proj_tif, filename=paste(paste(paste(path, task, ssp, scen, isos[which(iso[k]==isos)], sep = "/"), "pop_urb", ssp, scen, t, p,
  #                                               sep = "_"), "tif", sep = "."), format="GTiff", bylayer = T, overwrite=TRUE)
  #writeRaster(rur_proj_tif, filename=paste(paste(paste(path, task, ssp, scen, isos[which(iso[k]==isos)], sep = "/"), "pop_rur", ssp, scen, t, p,
  #                                               sep = "_"), "tif", sep = "."), format="GTiff", bylayer = T, overwrite=TRUE)
  
  list(stack(tot_proj_tif, urb_proj_tif, rur_proj_tif), err)

} #end of k loop

tifs <- results[[1]]
err <- results[[2]]

# make df of error results
err_df <- data.frame(col_iso=character(),
                     abs_sum <- numeric(),
                     abs_mn <- numeric(),
                     rtae <- numeric(),
                     wmape <- numeric(),
                     rmse <- numeric(),
                     rmse_p <- numeric(),
                     rsq <- numeric(),
                     stringsAsFactors=FALSE) 

for (i in 1:length(err)) {
  err_cntr <- err[[i]]
  col_iso <- rep(iso[i], nrow(err_cntr))
  err_cntr <- data.frame(col_iso, err_cntr)
  err_df <- rbind(err_df, err_cntr)
}
colnames(err_df) <-  c("iso", "abs_sum", "abs_mn", "rtae", "wmape", "rmse", "rmse_p", "rsq")


write.csv(err_df, paste0(paste(paste(path, task, "err_esti", sep = "/"), reg, "cntr", a, model_period, p, sep = "_"), ".csv"), row.names = F)


# ----------------------------- #
#### Step 5 - mosaic rasters ####
# ----------------------------- #

tot_proj_tif <- list()
urb_proj_tif <- list()
rur_proj_tif <- list()

for (i in 1:length(tifs)) {
  tot_proj_tif[[i]] <- tifs[[i]][[1:2]]
  urb_proj_tif[[i]] <- tifs[[i]][[3:4]]
  rur_proj_tif[[i]] <- tifs[[i]][[5:6]]
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

# write rasters as tif
writeRaster(tot_reg, filename=paste(paste(paste(path, task, reg, sep = "/"), "pop_tot", a, t, p, sep = "_"), "tif", 
                                    sep = "."), format="GTiff", bylayer = T, overwrite=TRUE)
writeRaster(urb_reg, filename=paste(paste(paste(path, task, reg, sep = "/"), "pop_urb", a, t, p, sep = "_"), "tif", 
                                    sep = "."), format="GTiff", bylayer = T, overwrite=TRUE)
writeRaster(rur_reg, filename=paste(paste(paste(path, task, reg, sep = "/"), "pop_rur", a, t, p, sep = "_"), "tif", 
                                    sep = "."), format="GTiff", bylayer = T, overwrite=TRUE)


# save entire workspace
save.image(paste(paste(paste(path, task, task, sep = "/"), reg, a, model_period, p, sep = "_"), "RData", sep = "."))

