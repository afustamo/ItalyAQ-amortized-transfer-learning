#Libraries
library(LatticeKrig)
library(imager)

setwd("C:/Users/anton/Desktop/Research/Italy_AQ_paper/Italy_AQ_AmortisedLatticeKrig")

# Load in Ale datasets
# cams_df <- readRDS("data/cams_we_em_2023.rds")
cams_df <- readRDS("data/cams_we_em_2023augmented.rds")
pred_df <- readRDS("data/pred_grid_df.rds")
stations_df <- readRDS("data/aq_eea_df.rds")
load("data/Station_registry_information.rda")
stations_coords <- Station_registry_information
rm(Station_registry_information)

# # Reorder cams_df because of lat issue
# cams_df <- cams_df[order(
#   cams_df$time,
#   cams_df$Latitude,        
#   cams_df$Longitude        
# ), ]

# cams_df_new <- readRDS("data/cams_we_em_2023augmented.rds")
# cams_df_new_crop <- cams_df_new[cams_df_new$time >= 19358 & cams_df_new$time <= 19722, ]
# 
# test.for.zero(cams_df$NO2, cams_df_new_crop$NO2)
# test.for.zero(cams_df$t2m, cams_df_new_crop$t2m)
# rm(cams_df, cams_df_new_crop)


rows_small <- 130
cols_small <- 130

gridList_small <- list( 
  x = seq( 6.0,18.9,by=.1), 
  y= seq( 35.0,47.9,by=.1)
)

sGrid_small<- make.surface.grid(gridList_small)

rows_big <- 259
cols_big <- 259

gridList_big <- list(
  x = seq( 6.0,18.9,by=.05), 
  y= seq( 35.0,47.9,by=.05)
)

sGrid_big<- make.surface.grid(gridList_big)

# Initial data visuals 
imagePlot(as.surface(gridList_small, cams_df$NO2[cams_df$time == 19358]), col=plasma(256), 
        main = "small"  ) 
imagePlot(as.surface(gridList_big, pred_df$NO2[pred_df$time == 19358]), col=plasma(256), 
          main = "big")

imagePlot(as.surface(gridList_big, pred_df$DEM[pred_df$time == 19358]), col=cividis(256), 
          main = "big DEM")

# should be 130 x 130 = 16900
length(cams_df$NO2[cams_df$time == 19358])
# should be 259 x 259 = 67081
length(pred_df$NO2[pred_df$time == 19358])

# DEM is missing bottom row (259)
sum(is.na(pred_df$DEM[pred_df$time == 19358]))


# ------------------------------------
# cropping CAMS down to 128x128 each 
# ------------------------------------

# Get unique lat/lon values that define the grid
lat_vals <- sort(unique(cams_df$Latitude))
lon_vals <- sort(unique(cams_df$Longitude))

nlat <- length(lat_vals)  # should be 130
nlon <- length(lon_vals)  # should be 130

# Assign row/col indices on the grid
cams_df$row_idx <- match(cams_df$Latitude,  lat_vals)  # 1 (north) ... 130 (south)
cams_df$col_idx <- match(cams_df$Longitude, lon_vals)  # 1 (west)  ... 130 (east)

cams_df <- subset(
  cams_df,
  row_idx >= 3 &       
    col_idx >= 2 &       
    col_idx <= 129       
)

# Optional: drop the helper indices if you don’t want them anymore
cams_df$row_idx <- NULL
cams_df$col_idx <- NULL

# should be 128 
sqrt(nrow(cams_df)/length(unique(cams_df$time)))

# should be 128 x 128 = 16384
length(cams_df$NO2[cams_df$time == 19358])

# adjust gridlist
gridList_small <- list( 
  x = seq( 6.1,18.8,by=.1), 
  y= seq( 35.2,47.9,by=.1)
)

sGrid_small<- make.surface.grid(gridList_small)

imagePlot(as.surface(gridList_small, cams_df$NO2[cams_df$time == 19358]), col=plasma(256), 
          main = "small cropped"  ) 


# ------------------------------------
# cropping PRED down to 255x255 each 
# ------------------------------------
pred_df$Longitude[pred_df$time == 19358]
pred_df$Latitude[pred_df$time == 19358]

length(pred_df$Longitude[pred_df$time == 19358])

lat_seq <- seq( 35.0,47.9,by=.05)
lon_seq <- seq( 6.0,18.9,by=.05)

grid_temp <- expand.grid(Longitude = lon_seq,
                         Latitude = lat_seq)
grid <- grid_temp[order(grid_temp$Latitude, grid_temp$Longitude), ]

times <- unique(pred_df$time)

for (tt in times) {
  idx <- which(pred_df$time == tt)
  
  stopifnot(length(idx) == 259*259)
  
  # Assign correct lat/lon
  pred_df$Latitude[idx]  <- grid$Latitude
  pred_df$Longitude[idx] <- grid$Longitude
}


pred_df$Longitude[pred_df$time == 19358]
pred_df$Latitude[pred_df$time == 19358]



# Get unique lat/lon values that define the grid
lat_vals <- sort(unique(pred_df$Latitude))
lon_vals <- sort(unique(pred_df$Longitude))

nlat <- length(lat_vals) 
nlon <- length(lon_vals)  

# Assign row/col indices on the grid
pred_df$row_idx <- match(pred_df$Latitude,  lat_vals)  
pred_df$col_idx <- match(pred_df$Longitude, lon_vals)  


pred_df <- subset(
  pred_df,
  row_idx >= 5 &       
    col_idx >= 3 &       
    col_idx <= 257       
)

# Optional: drop the helper indices if you don’t want them anymore
pred_df$row_idx <- NULL
pred_df$col_idx <- NULL

# should be 255
sqrt(nrow(pred_df)/365)

# should be 255 x 255 = 65025
length(pred_df$NO2[pred_df$time == 19358])

# adjust gridlist
gridList_big <- list( 
  x = seq( 6.1,18.8,by=.05), 
  y= seq( 35.2,47.9,by=.05)
)

par(mfrow = c(1,2))
imagePlot(as.surface(gridList_big, pred_df$NO2[pred_df$time == 19358]), col=plasma(256), 
          main = "big cropped"  ) 

imagePlot(as.surface(gridList_small, cams_df$NO2[cams_df$time == 19358]), col=plasma(256), 
          main = "small cropped"  ) 
par(mfrow = c(1,1))

# ----------------------------------------
# downscaling DEM and adding it to cams_df
# ----------------------------------------
# 
# # practice run on one 
# elev_map <- pred_df$DEM[pred_df$time == 19358]
# elev_mat <- matrix(elev_map, nrow = 255, ncol = 255)
# 
# img <- as.cimg(elev_mat)
# 
# elev_small <- resize(img, size_x = 128, size_y = 128, interpolation_type = 2)
# elev_small_mat <- as.matrix(elev_small)
# 
# par(mfrow = c(1,2))
# imagePlot(elev_small_mat, col = cividis(256),
#           main = "Downsampled DEM (128 x 128)")
# imagePlot(elev_mat, col = cividis(256),
#           main = "Original DEM (255 x 255)")
# par(mfrow = c(1,1))
# 
# 
# # now let's do it to all of them 
# cams_df$DEM <- NA
# 
# # Grid sizes
# nx_old <- 255
# ny_old <- 255
# nx_new <- 128
# ny_new <- 128
# 
# # All timesteps
# times <- sort(unique(pred_df$time))
# nt    <- length(times)
# 
# # Pre-allocate: one vector per time, each of length 128*128
# dem_list <- vector("list", nt)
# 
# for (k in seq_along(times)) {
#   tt <- times[k]
#   
#   # 1) Extract DEM for this time
#   elev_map <- pred_df$DEM[pred_df$time == tt]
#   
#   # 2) Make 255 x 255 matrix
#   elev_mat <- matrix(elev_map, nrow = nx_old, ncol = ny_old)
#   
#   # 3) Convert to cimg and resize to 128 x 128
#   img <- as.cimg(elev_mat)
#   elev_small <- resize(img,
#                        size_x = nx_new,
#                        size_y = ny_new,
#                        interpolation_type = 2)  # your chosen method
#   
#   # 4) Back to matrix → vector (row-major)
#   elev_small_mat <- as.matrix(elev_small)
#   dem_list[[k]]  <- as.vector(elev_small_mat)
# }
# 
# # Combine all times into one long vector: 128*128*nt
# dem_128_vec <- unlist(dem_list)
# 
# # Sanity check: cams_df must match this length
# stopifnot(nrow(cams_df) == length(dem_128_vec))
# 
# # 5) Append to cams_df
# cams_df$DEM <- dem_128_vec
# rm(dem_128_vec)
# rm(dem_list)
# rm(grid_temp)

sample_time <- 19378

par(mfrow = c(1,2))
imagePlot(as.surface(gridList_big, pred_df$DEM[pred_df$time == sample_time]), 
          col=cividis(256), 
          main = "Original DEM (255 x 255)"  )
imagePlot(as.surface(gridList_small, cams_df$DEM[cams_df$time == sample_time]), 
          col=cividis(256), 
          main = "Downscaled DEM (128 x 128)"  )
par(mfrow = c(1,1))


# remove predicted lm/arx1 columns
# cams_df <- cams_df[, c(1:14, 21)]
pred_df <- pred_df[,c(1:14, 21,22)]


# ------------------------------------
# adding lat lon to stations 
# ------------------------------------

stations_df <- merge(
  stations_df,
  stations_coords[, c("AirQualityStation", "Latitude", "Longitude")],
  by = "AirQualityStation",
  all.x = TRUE
)

stations_df <- stations_df[,c(1:14, 21:24)]
rm(stations_coords)

# save final data
saveRDS(cams_df, file = "data/cams_df_2023padded.rds")
saveRDS(pred_df, file = "data/pred_df_2023.rds")
saveRDS(stations_df, file = "data/eea_df_2023.rds")

