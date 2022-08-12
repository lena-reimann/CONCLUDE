###############################
## Produce distance matrixes ##
###############################
# by Lena Reimann
# Jan 05, 2021

# Goal: produce distance matrix for each cntr (on HPC cluster)
# comments: # 1) use as input for beta optimization, alpha calculation, and population projections


rm(list=ls())

# load packages 
#.libPaths("~/r_libs")

library(raster, lib.loc = "~/r_libs")  
library(ggplot2, lib.loc = "~/r_libs") 
library(geosphere, lib.loc = "~/r_libs") 
library(scales, lib.loc = "~/r_libs") 
library(Matrix, lib.loc = "~/r_libs")
library(data.table, lib.loc = "~/r_libs")

#for parallel processing
library("parallel", lib.loc = "~/r_libs") 
library("foreach", lib.loc = "~/r_libs") 
library("doParallel", lib.loc = "~/r_libs") 
library("igraph", lib.loc = "~/r_libs") 
library(doSNOW, lib.loc = "~/r_libs") 

chunks <- 4  # define the number of cores to use in parallel (8 for cluster; 4 for PC)

#-------------------------------------------#
# Step 1 - Define basic information of task #
#-------------------------------------------#
task            <- "dist_m"
reg             <- "N" # define region
id              <- c(705) # define cntr ids to be used
iso             <- c("SVN") # define cntr ISOs to be used
gw              <- 20      # gravity window used
earth_r         <- 6378.137  # Earth's radius in km


#----------------------------------------#
# Step 2 - Load & preprocess data needed #
#----------------------------------------#

# define path where data are located
path <- "G:/PhD/Fulbright/Demographic_model/model_extension"
setwd(paste(path, sep = "/")) 
getwd()


# (1) define gpw ids
gpw_ids <- read.csv("boundary_data/gpw4_10/gpw_id_med_iso25.csv")
isos <- gpw_ids$iso
ids <- gpw_ids$gpw_id


# (2) load mask data 
#per cntr
mask_cntr <- list()
#load masks per country for each t
for (i in seq(ids)) {
  name <- paste0(paste(paste0("spatial_mask/sp_mask_cntr25_ssp3/sp_mask_ssp3_id", ids[i]), "2000", sep = "_"), ".tif")  
  mask_cntr[[i]] <- raster(name)
}


# (3) load distance raster for reducing the number of grid cells of large countries
buf <- raster("pop_data/pop_2015/pop15_R2018_eucl.tif") # distance raster produced in ArcGIS
buf[buf > 30000] <- NA
buf[!is.na(buf)] <- 1


# preprocessing: produce lat/lon coordinates per cntr needed for running the functions
cntr_LL <- list()

for (k in 1:length(ids)) {
  ## mask ##
  # extract mask data needed
  mask00 <- mask_cntr[[which(ids[k]==ids)]]
  
  ## for very large cntr with limited area available for settlement (i.e. dza, lby) ## 
  #reduce number of cells to be processed based on a buffer layer
  if (ids[k] == 12 || ids[k] == 434 || ids[k] == 818) { # if dza or lby to be projected
    buf30 <- resample(buf, mask00, method = "ngb")
    mask <- mask00
    
    # if dza reduce the buffer to 20 km to reduce processing time
    if (ids[k] == 12) {
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
  
  cntr_LL[[k]] <- mask[,c(1,2)] # get lon/lat based on the spatial mask
  
}


#---------------------------#
# Step 3 - define functions #
#---------------------------#
# Create a sparse distance matrix, taking into account the gravity window whose application
# will result in a large number of zero distances.

# (1) compute distances (to be used in 2.)
compute_distances <- function(country_cntr_LL, earth_r, gravity_window, row_start, row_end) {
  print(paste("Processing: start row =", row_start))
  print(paste("Processing: end row =", row_end))
  num_cells <- nrow(country_cntr_LL)
  # Holds the column indices of the non-zero elements of the matrix
  nz_col_indices <- list()
  # Holds the row indices of the non-zero elements of the matrix
  nz_row_indices <- list()
  # Holds the actual non zero values
  nz_data <- list()
  
  list_index = 1
  msg_cnt <- floor((row_end - row_start + 1) * 0.05)
  iter_cnt <- 1
  t_start <- Sys.time()
  for ( i in row_start:row_end ) {
    if ( msg_cnt > 0 && iter_cnt %% msg_cnt == 0) {
      #print(paste("Processed :", iter_cnt / (row_end - row_start + 1)))
      t_stop <- Sys.time()
      #print(t_stop - t_start)
      t_start <- t_stop
    }
    i_point = country_cntr_LL[i, ]
    # Get the distances of point i with the rest (taking into account the symmetry of the distances).
    dist_pairs_irow <- distHaversine(i_point, country_cntr_LL[i:num_cells,], earth_r)
    non_zero_j_indices <- which(dist_pairs_irow <= gravity_window & dist_pairs_irow > 0 )
    
    # If there are diastances within the gravity window update the non zero column, rows and data
    # vectors accordingly.
    if ( length(non_zero_j_indices) > 0 ) {
      nz_col_indices[[list_index]] <- (non_zero_j_indices + (i - 1))
      nz_row_indices[[list_index]] <- c(rep(i, length(non_zero_j_indices)))
      nz_data[[list_index]] <- dist_pairs_irow[non_zero_j_indices]
      list_index <- list_index + 1
    }
    iter_cnt <- iter_cnt + 1
  }
  
  temp_list <- list("nz_rows" = as.vector(unlist(nz_row_indices)), 
                    "nz_cols" =as.vector(unlist(nz_col_indices)), 
                    "nz_data" = as.vector(unlist(nz_data)))
  return(as.data.frame(temp_list))
}

# (2) calc distance matrix
get_distance_matrix_parallel <- function(country_cntr_LL, earth_r, gravity_window, num_chunks) {
  num_rows <- nrow(country_cntr_LL)
  chunk_size = floor(num_rows / num_chunks)
  chunk_rem = num_rows %% num_chunks
  
  sparse_matrix_info <- foreach(chunk = 1:num_chunks, 
                                .combine = "rbind", 
                                .packages = c("raster", "igraph", "geosphere"),
                                .export = 'compute_distances',
                                .verbose = TRUE) %dopar% {
                                  row_start <- (chunk - 1) * chunk_size + 1
                                  row_end <- chunk * chunk_size
                                  chunk_matrix <- compute_distances(country_cntr_LL, earth_r, gravity_window, row_start, row_end)
                                }
  
  if ( chunk_rem > 0) {
    row_start <- num_chunks * chunk_size
    row_end <- num_rows
    rem_chunk_matrix <- compute_distances(country_cntr_LL, earth_r, gravity_window, row_start, row_end)
    sparse_matrix_info <- rbind(sparse_matrix_info, rem_chunk_matrix)
  }
  
  return(sparseMatrix(i = sparse_matrix_info$nz_rows,
                      j = sparse_matrix_info$nz_cols,
                      x = sparse_matrix_info$nz_data,
                      dims = c(num_rows, num_rows), symmetric = TRUE))
}

# (3) function for writing compressed file of MM matrix
# @param x A sparse matrix from the Matrix package.
# @param file A filename that ends in ".gz".
writeMMgz <- function(x, file) {
  mtype <- "real"
  if (is(x, "ngCMatrix")) {
    mtype <- "integer"
  }
  writeLines(
    c(
      sprintf("%%%%MatrixMarket matrix coordinate %s symmetric", mtype),
      sprintf("%s %s %s", x@Dim[1], x@Dim[2], length(x@x))
    ),
    gzfile(file)
  )
  data.table::fwrite(
    x = summary(x),
    file = file,
    append = TRUE,
    sep = " ",
    row.names = FALSE,
    col.names = FALSE
  )
}


#----------------------------------------------------------#
# Step 4 - Run functions for producing the distance matrix #
#----------------------------------------------------------#

# define path for saving the results (in case it deviates from the one defined above)
path <- "G:/PhD/Fulbright/Demographic_model/model_extension/Theo" 

cl <- makeCluster(chunks)
registerDoParallel(cl, cores = chunks)
t_start_all <- Sys.time()

for (i in 1:length(id)) {
dist_sp <- get_distance_matrix_parallel(cntr_LL[[which(id[i]==ids)]], earth_r, gw, chunks)

# save
name <- paste(paste(path, task, isos[which(id[i]==ids)], sep = "/"), "mtx", sep = ".")
writeMM(dist_sp, name) #uncompressed

name <- paste(paste(path, task, isos[which(id[i]==ids)], sep = "/"), "mtx.gz", sep = ".")
writeMMgz(dist_sp, name) # compressed

}
t_stop_all <- Sys.time()
print(t_stop_all - t_start_all)

stopCluster(cl) 

















