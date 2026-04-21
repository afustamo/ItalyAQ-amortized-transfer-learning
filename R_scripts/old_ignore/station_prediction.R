#Libraries
library(rhdf5)
library(LatticeKrig)
library(sf)

#import hacked LK files for additional kappa2 estimation
setwd("C:/Users/anton/Desktop/Research/Italy_AQ_paper")
source("Italy_AQ_AmortisedLatticeKrig/R_scripts/hacked_LatticeKrig.R")

set.seed(7)

# create grid 
rows_small <- length(seq( 6.1,18.8,by=.1))

gridList_small <- list( 
  x = seq( 6.1,18.8,by=.1), 
  y= seq( 35.2,47.9,by=.1)
)
sGrid_small<- make.surface.grid(gridList_small)

# sDomain <- cbind( range(gridList_small$x), range(gridList_small$y))

# need this predict function for non-stationary, anisotropic awght
predict.multivariateSurfaceGrid<- function(object,x){
  dimZ<- dim( object$z)
  L<- dimZ[3]
  out<- matrix( NA, nrow= nrow(x), ncol=L)
  for (  l in 1:L){
    out[,l]<- interp.surface( 
      list( x=object$x,y=object$y, z=object$z[,,l]) , x)
  }
  return( out)
}

# ------------------------------------
# Function to create an LKinfo object
# from STUN parameters
# ------------------------------------

# simply specify a file path, dataset name and
# a day (1-365) in 2023
# different datasets of parameters with different normalizations
# which were processed with different replicate counts

# dataset name options are: 
# "lm_rscale_1rep_output"
# "arx1_rscale_1rep_output"
# "arx1_surround_1rep_output" 
# "arx1_rscale_30rep_output"
# "arx1_surround_30rep_output"

# different datasets of parameters with different normalizations
# which were processed with different replicate counts

make_nonstat_LKinfo <- function(
    file_path, 
    dataset_name,
    day, 
    gridlist,
    normalize = TRUE, 
    NC = 128, 
    nlevel = 1, 
    NC.buffer = 0, 
    sanity_plotting = FALSE, 
    sanity_sim = FALSE
){
  params <- h5read(file_path, dataset_name)
  H5close()
  # for plotting and image purposes, need to flip upside down
  params <- params[,128:1,,day]
  
  # recover kappa
  kappa2 <- exp(params[,,1])
  # just in case we need awght
  awght <- kappa2 + 4
  # theta needs to be transformed like this 
  theta <- params[,,2] + pi/2
  rho <- params[,,3]
  
  if (sanity_plotting){
    par(mfrow = c(1,3))
    imagePlot(as.surface(gridlist, kappa2), main="kappa2", col = viridis(256))
    world(add=TRUE, col = "white", lwd = 1)
    imagePlot(as.surface(gridlist, theta-pi/2), main="theta (adj)", col = viridis(256))
    world(add=TRUE, col = "white", lwd = 1)
    imagePlot(as.surface(gridlist, rho), main="rho", col = viridis(256))
    world(add=TRUE, col = "white", lwd = 1)
    par(mfrow = c(1,1))
  }
  
  # need these for encoding into LK
  rhox <- sqrt(rho)
  rhoy <- 1/rhox
  
  # create H tensor out of params
  H11 <- ( rhox^2 * (cos(theta))^2) + ( rhoy^2 * (sin(theta))^2 ) 
  H12 <- (rhoy^2 - rhox^2)*(sin(theta)*cos(theta))
  H21 <- H12 
  H22 <- (rhox^2 * (sin(theta))^2) + (rhoy^2 * (cos(theta))^2)
  
  rows <- length(gridlist$x)
  
  # fill the high dimensional stencil (9 fields)
  stencil_tensor <- array( NA, c( rows,rows,9))
  stencil_tensor[,,1] <- 0.5 * H12
  stencil_tensor[,,2] <- -H22
  stencil_tensor[,,3] <- -0.5 * H12
  stencil_tensor[,,4] <- -H11
  stencil_tensor[,,5] <- kappa2 + 2 * H11 + 2 * H22
  stencil_tensor[,,6] <- -H11
  stencil_tensor[,,7] <- -0.5 * H12
  stencil_tensor[,,8] <- -H22
  stencil_tensor[,,9] <- 0.5 * H12
  
  # next, we put everything into awght obj of a particular class
  awght_obj <- list( x= gridlist$x,  y= gridlist$y, z=stencil_tensor )
  class( awght_obj)<- "multivariateSurfaceGrid"
  
  sGrid <- make.surface.grid(gridlist)
  
  LKinfo <- LKrigSetup(
    # dont change grid
    sGrid, 
    # you change awght indirectly with day and file selection
    a.wghtObject =  awght_obj,
    # the rest of these you change directly in function calls 
    NC = NC, 
    nlevel = nlevel,
    normalize=normalize,
    NC.buffer = NC.buffer
  )
  
  if(sanity_sim){
    test <- LKrig.sim(
      sGrid,
      LKinfo=LKinfo,
      M = 1
    )
    
    if (sanity_plotting){
      imagePlot(as.surface(gridlist, test), col = turbo(256))
      world(add=TRUE, col = "black", lwd = 1)
    }
  }
  
  return(LKinfo)
}

# choose a day for this analysis (24, 77, 178, 256, 319)
chosen_day <- 256

# example usage 
LKinfo_ns <- make_nonstat_LKinfo(
  file_path = "Italy_AQ_AmortisedLatticeKrig/data/STUN_param_df_2023.h5",
  dataset_name = "arx1_surround_30rep_output",
  day = chosen_day,
  gridlist = gridList_small,
  normalize = TRUE,
  sanity_plotting = TRUE, 
  sanity_sim = TRUE, 
)

# creating our stationary model without the func
LKinfo_stat <- LKrigSetup(
  sGrid_small,
  NC = 128, 
  nlevel = 1,
  normalize=TRUE,
  NC.buffer = 0, 
  a.wght = 4.1
)



# ------------------------------------
# 1. Is STUN a good estimator? 
# Synthetic experiment with CTM 
# values at station locations
# NEED TO REDO THIS ENTIRELY
# ------------------------------------
stations_df <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/eea_df_2023.rds")

# counting non-NAs. remove unname if u want dates
# non_na_counts <- unname(
#   tapply(
#     stations_df$EEA_NO2, stations_df$time, function(x) sum(!is.na(x))
#   )
# )
# non_na_counts

# subset down to our chosen day 
stations_df_day <- stations_df[
  stations_df$time == (stations_df$time[1] + chosen_day -1),
]
rm(stations_df)

cams_df <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/cams_df_2023.rds")
cams_df_day <- cams_df[
  cams_df$time == (cams_df$time[1] + chosen_day -1),
]
rm(cams_df)

cams_df_day[,"ssr"] <- cams_df_day[,"ssr"] / 10^6
stations_df_day[,"ssr"] <- stations_df_day[,"ssr"] / 10^6


stations_sf_day <- st_as_sf(
  stations_df_day, 
  coords = c("Longitude", "Latitude"), 
  crs = 4326
)
cams_sf_day <- st_as_sf(
  cams_df_day, 
  coords = c("Longitude", "Latitude"), 
  crs = 4326
)

nearest_idx <- st_nearest_feature(stations_sf_day, cams_sf_day)
train_idx <- unique(nearest_idx)

# train_idx <- sample(1:nrow(cams_df_day), size = floor(0.005*nrow(cams_df_day)), replace = FALSE)

train_df <- cams_df_day[train_idx, ]
test_df <- cams_df_day[-train_idx, ]

s <- train_df[, c("Longitude", "Latitude")]
y <- train_df$NO2
Z <- as.matrix(
  train_df[, c("rh", "ssr", "t2m",
                      "windspeed", "sl_blh", "EM_NO2", "lag_cams_no2")]
)

par(mfrow = c(1,2))
bubblePlot(s, y, main = "Train", col = turbo(256))

bubblePlot(test_df[, c("Longitude", "Latitude")], 
           test_df$NO2, main = "Test", col = turbo(256))
par(mfrow = c(1,1))

linear_model <- lm(
  NO2 ~ Longitude + Latitude + rh + ssr + t2m + windspeed + sl_blh + EM_NO2 + lag_cams_no2, 
  data = train_df
)

Lmodel_stat <- LatticeKrig(
  s,
  y,
  Z = Z,
  LKinfo = LKinfo_stat,
  findAwght = TRUE
)

Lmodel_ns <- LatticeKrig(
  x = s,
  y = y,
  Z = Z,
  LKinfo = LKinfo_ns
)

linear_model
Lmodel_stat
Lmodel_ns

summary(abs(Lmodel_stat$residuals))
summary(abs(Lmodel_ns$residuals))
summary(abs(linear_model$residuals))


# prediction on test set: 
s_test <- test_df[, c("Longitude", "Latitude")]
y_test <- test_df$NO2
Z_test <- as.matrix(
  test_df[, c("rh", "ssr", "t2m",
                  "windspeed", "sl_blh", "EM_NO2", "lag_cams_no2")]
)

pred_linear <- predict(linear_model, newdata = test_df)
pred_stat <- predict(Lmodel_stat, s_test, Z = Z_test)
pred_ns <- predict(Lmodel_ns, s_test, Z = Z_test)

# for separate OLS
# pred_linear <- predict(linear_model, newdata = test_df)
# pred_stat <- predict(Lmodel_stat, s_test) + pred_linear
# pred_ns <- predict(Lmodel_ns, s_test) + pred_linear

sqrt(mean((pred_linear - y_test)^2))
sqrt(mean((pred_stat - y_test)^2))
sqrt(mean((pred_ns - y_test)^2))

summary(abs(pred_stat - y_test))
summary(abs(pred_ns - y_test))




# # now let's predict
Z_map <- as.matrix(
  cams_df_day[, c("rh", "ssr", "t2m",
                  "windspeed", "sl_blh", "EM_NO2", "lag_cams_no2")]
)

cams_grid <- make.surface.grid(gridList_small)

surf_pred_stat <- predict(Lmodel_stat, cams_grid, Z = Z_map)
surf_pred_stat_mu <- predict(Lmodel_stat, cams_grid, Z = Z_map, just.fixed = TRUE)
surf_pred_stat_gp <- surf_pred_stat - surf_pred_stat_mu

surf_pred_ns <- predict(Lmodel_ns, cams_grid, Z = Z_map)
surf_pred_ns_mu <- predict(Lmodel_ns, cams_grid, Z = Z_map, just.fixed = TRUE)
surf_pred_ns_gp <- surf_pred_ns - surf_pred_ns_mu

sqrt(mean((surf_pred_stat - cams_df_day$NO2)^2))
sqrt(mean((surf_pred_ns - cams_df_day$NO2)^2))





# for separate OLS 
surf_pred_linear <- predict(linear_model, newdata = cams_df_day)
sqrt(mean((surf_pred_linear - cams_df_day$NO2)^2))



par(mfrow = c(2,3))
par(mar = c(2,2,2,2))
plotmax <- max(
  c(surf_pred_stat, surf_pred_ns, cams_df_day$NO2)
)
plotmin <- min(
  c(surf_pred_stat, surf_pred_ns, cams_df_day$NO2)
)

imagePlot(as.surface(gridList_small, surf_pred_stat), 
          main="Stationary LK prediction", col = turbo(256), 
          zlim = c(plotmin, plotmax))
# world(add=TRUE, col = "white", lwd = 1)
imagePlot(as.surface(gridList_small, surf_pred_ns), 
          main="Non-stationary LK prediction", col = turbo(256), 
          zlim = c(plotmin, plotmax))
# world(add=TRUE, col = "white", lwd = 1)

imagePlot(as.surface(gridList_small, cams_df_day$NO2), 
          main="CAMS NO2", col = turbo(256), 
          zlim = c(plotmin, plotmax))


res_max <- max(
  surf_pred_stat - cams_df_day$NO2, 
  surf_pred_ns - cams_df_day$NO2
)
res_min <- min(
  surf_pred_stat - cams_df_day$NO2, 
  surf_pred_ns - cams_df_day$NO2
)
imagePlot(as.surface(gridList_small, surf_pred_stat - cams_df_day$NO2), 
          main="Residuals: Stationary", col = turbo(256), zlim = c(res_min, res_max))
# world(add=TRUE, col = "white", lwd = 1)
imagePlot(as.surface(gridList_small, surf_pred_ns - cams_df_day$NO2), 
          main="Residuals: Non-stationary", col = turbo(256), zlim = c(res_min, res_max))
# world(add=TRUE, col = "white", lwd = 1)

imagePlot(as.surface(gridList_small,surf_pred_ns - surf_pred_stat), 
          main="Difference (Nonstat - stat)", col = turbo(256))


par(mfrow = c(2,3))

imagePlot(as.surface(gridList_small, surf_pred_stat_mu), 
          main="Stationary LK mean", col = turbo(256))
# world(add=TRUE, col = "white", lwd = 1)
imagePlot(as.surface(gridList_small, surf_pred_ns_mu), 
          main="Non-stationary LK mean", col = turbo(256))
# world(add=TRUE, col = "white", lwd = 1)

plotmin <- min(
  cams_df_day$arx1_resid - surf_pred_stat_mu, 
  surf_pred_stat_gp, surf_pred_ns_gp
)
plotmax <- max(
  cams_df_day$arx1_resid - surf_pred_stat_mu, 
  surf_pred_stat_gp, surf_pred_ns_gp
)

imagePlot(as.surface(gridList_small, cams_df_day$arx1_resid - surf_pred_stat_mu), 
          main="Fixed Effect Residual", col = turbo(256), zlim = c(plotmin, plotmax))

imagePlot(as.surface(gridList_small, surf_pred_stat_gp), 
          main="Stationary LK GP", col = turbo(256), zlim = c(plotmin, plotmax))
# world(add=TRUE, col = "white", lwd = 1)
imagePlot(as.surface(gridList_small, surf_pred_ns_gp), 
          main="Non-stationary LK GP", col = turbo(256), zlim = c(plotmin, plotmax))
# world(add=TRUE, col = "white", lwd = 1)

imagePlot(as.surface(gridList_small,surf_pred_ns_gp - surf_pred_stat_gp), 
          main="GP Difference (NS - Stat)", col = turbo(256))

par(mfrow = c(1,1))
