##########             IMPORTANT              ###########
#-------------------------------------------------------#
# Make sure to set your working directory 
# "To Source File Location". 
# Also, make sure you have downloaded all required 
# packages from required_packages.R. 
# Now you can run the script. 
#-------------------------------------------------------#

#Libraries
library(rhdf5)
library(LatticeKrig)
library(tictoc)
library(here)

#---------------------------------
# Fitting ARX1 to data 
#---------------------------------

cams_df <- readRDS(here("data", "cams_df_2023padded.rds"))
# first day does not have lagged NO2 (19336)
cams_df <- cams_df[cams_df$time >= 19337, ] 

#split by day 
cams_by_day <- split(cams_df, cams_df$time)

# fit an ARX(1) model for each day
arx1_by_day <- lapply(
  cams_by_day,
  function(df_day) {
    lm(
      NO2 ~ Longitude + Latitude + rh + ssr + t2m + windspeed + sl_blh +
        EM_NO2 + DEM + lag_cams_no2, data = df_day)
  }
)

rm(cams_by_day)

# add results to cams df 
cams_df$arx1_predict <- NA
cams_df$arx1_resid <- NA

#coefs to compare to eea model
arx1_coefs <- data.frame(
  Intercept = NA, 
  Longitude = NA,  
  Latitude  = NA, 
  rh  = NA, 
  ssr = NA, 
  t2m = NA, 
  windspeed = NA, 
  sl_blh = NA, 
  EM_NO2 = NA, 
  DEM = NA, 
  lag_cams_no2 = NA
)


i <- 1
for (tt in names(arx1_by_day)) {
  
  # convert the loop key to Date
  tt_date <- as.Date(tt)
  
  # match the df (make sure cams_df$time is Date too)
  idx <- as.Date(cams_df$time) == tt_date
  
  cams_df$arx1_predict[idx] <- arx1_by_day[[i]]$fitted.values
  cams_df$arx1_resid[idx]   <- arx1_by_day[[i]]$residuals
  
  arx1_coefs[i, ] <- arx1_by_day[[i]]$coefficients
  i <- i + 1
}



#---------------------------------
# Checking results with plotting
#---------------------------------

rows <- 128
cols <- 128

gridList <- list( 
  x = seq( 6.1,18.8,by=.1), 
  y= seq( 35.2,47.9,by=.1)
)
sGrid<- make.surface.grid(gridList)


# choose a time
time_idx <- 19358
time_idx <- 19460


# first, predictions
par(mfrow = c(1,2))
pred_max <- max(
  cams_df$NO2[cams_df$time == time_idx], 
  cams_df$arx1_predict[cams_df$time == time_idx]
)
pred_min <- min(
  cams_df$NO2[cams_df$time == time_idx], 
  cams_df$arx1_predict[cams_df$time == time_idx]
)
imagePlot(as.surface(gridList, cams_df$NO2[cams_df$time == time_idx]), col=plasma(256), 
          main = "Original",  zlim = c(pred_min, pred_max)) 

imagePlot(as.surface(gridList, cams_df$arx1_predict[cams_df$time == time_idx]), col=plasma(256), 
          main = "arx1 predict", zlim = c(pred_min, pred_max)) 
par(mfrow = c(1,1))


# now some residuals
par(mfrow = c(1,1))
res_max <- max(
  cams_df$arx1_resid[cams_df$time == time_idx]
)
res_min <- min(
  cams_df$arx1_resid[cams_df$time == time_idx]
)

imagePlot(as.surface(gridList, cams_df$arx1_resid[cams_df$time == time_idx]), col=turbo(256), 
          main = "arx1 resid")#, zlim = c(res_min, res_max)) 
par(mfrow = c(1,1))

# final save of data
saveRDS(cams_df, file = here("data", "cams_df_2023padded.rds"))
# 21 extra days in the beginning
# 15 extra days in the end

#---------------------------------
# saving residuals as hdf5 files 
#---------------------------------

#filepath to where the h5 file will be saved 
h5_file_path <- here("data", "NO2_resid_df_2023padded.h5")

# create the file
h5createFile(h5_file_path)

rows <- sqrt(nrow(cams_df)/length(unique(cams_df$time)))
cols <- rows

arx1_resid_forh5 <- array(
  data = cams_df$arx1_resid, 
  dim = c(rows, cols, length(unique(cams_df$time)))
)

h5write(obj = arx1_resid_forh5, file = h5_file_path, "arx1_resid") 

h5ls(h5_file_path)

H5close()

