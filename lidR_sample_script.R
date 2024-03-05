#### LIDR CHM PROCESSING SAMPLE ####

library(lidR) ; library(dplyr) ; library(terra) ; library(rgdal) ; library(raster) ; library(future)

# Read 50x50m .laz tiles as LASCatalog
# processing will occur by tile 
CTG = readLAScatalog(folder = "path_to_laz_processing_tiles")

##### GROUND CLASSIFICATION ####
# buffer each tile by 5m to reduce edge effects
opt_chunk_buffer(CTG) = 5
opt_laz_compression(CTG) = FALSE
# taking a random 10% of the point cloud speeds things up and reduces outliers 
opt_filter(CTG) = "-keep_random_fraction 0.1" 
# save classified pointcloud tiles to GROUND with corresponding tile ID 
opt_output_files(CTG) = paste0("\\output\\GROUND\\{*}_GROUND")  
opt_progress(CTG) = TRUE

# classify ground and canopy from the filtered point cloud using cloth simulation function
# set up and execute task in parallel (6 cores) 
plan(multisession, workers = 6L)

GROUND = classify_ground(CTG, algorithm = csf(
  # in case there are steep slopes 
  sloop_smooth = TRUE,
  class_threshold = 0.25,
  # larger resolution = coarser DTM
  cloth_resolution = 0.5,
  # rigidness of cloth - higher = more rigid
  rigidness = 3L,
  iterations = 500L,
  time_step = 1)) 

# always stop the parallel processing!! 
future:::ClusterRegistry("stop")

#### CREATE DEM FROM GROUND CLASSIFIED POINTS ####

# Read in the ground classified tiles
GROUND = readLAScatalog(folder = paste0(dir, "\\output\\GROUND\\"))
opt_chunk_buffer(GROUND) = 1
# drop remaining noisy points below the ground
opt_filter(CTG) = "-drop_z_below -0.25"
opt_laz_compression(CTG) = TRUE
opt_output_files(GROUND) = ""
opt_progress(GROUND) = TRUE

# Create a DTM from the classified ground points using 
# a k-nearest neighbour (KNN) approach with an inverse-distance weighting (IDW)
plan(multisession, workers = 6L)
DTM = grid_terrain(GROUND, res = 1, knnidw(), full_raster = TRUE)
future:::ClusterRegistry("stop")

# save the DTM 
writeRaster(DTM, "DTM_filepath.tif", overwrite = TRUE) 

#### NORMALIZE THE CANOPY POINTS TO THE DEM ####

# read the tiles again
CTG = readLAScatalog(folder = paste0(dir, "\\", site,"_laz\\"))
opt_chunk_buffer(CTG) = 1
opt_laz_compression(CTG) = TRUE
# regularize point cloud with a small 1cm voxel
opt_filter(CTG) = "-thin_with_voxel 0.01"
# write tiles to NORM folder 
opt_output_files(CTG) = paste0(dir, "\\output\\NORM\\{*}_NORM")
opt_progress(CTG) = TRUE

# normalize the point cloud to the DTM
plan(multisession, workers = 6L)
NORM = normalize_height(CTG, DTM, na.rm = TRUE)
future:::ClusterRegistry("stop")

# remove normalized points below the ground
# now read them in again and filter out points below -0.25m and above 30m from normalized surface
NORM = readLAScatalog(paste0(dir, "\\output\\NORM\\"))
# no buffer, just filtering
opt_chunk_buffer(NORM) = 0
opt_laz_compression(NORM) = TRUE
opt_filter(NORM) = "-drop_z_below -.1 -drop_z_above 30"
opt_output_files(NORM) = paste0(dir, "\\output\\NORM\\{*}")
opt_progress(NORM) = TRUE

# overwrite over the NORM files with the cleaned point clouds
plan(multisession, workers = 6L)
NORM_clean = catalog_retile(NORM)
future:::ClusterRegistry("stop")

#### CREATE A RASTERIZED CHM ####
# read in cleaned norm cloud
NORM = readLAScatalog(paste0(dir, "\\output\\NORM\\"))
opt_chunk_buffer(NORM) = 5
opt_laz_compression(NORM) = TRUE
opt_output_files(NORM) = ""
opt_progress(NORM) = TRUE

# create a 10cm CHM raster using maximum point value in pixel
plan(multisession, workers = 6L)
CHM_max = pixel_metrics(NORM,
                        res = 0.1,
                        func = ~max(Z))
future:::ClusterRegistry("stop")