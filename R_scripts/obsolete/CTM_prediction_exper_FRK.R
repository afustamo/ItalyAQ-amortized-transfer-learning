# Libraries
library(rhdf5)
library(LatticeKrig)
library(sf)

# import hacked LK files for additional kappa2 estimation
setwd("C:/Users/anton/Desktop/Research/Italy_AQ_paper")

set.seed(7)

# create grid 
rows_small <- length(seq(6.1, 18.8, by = .1))

gridList_small <- list(
  x = seq(6.1, 18.8, by = .1),
  y = seq(35.2, 47.9, by = .1)
)
sGrid_small <- make.surface.grid(gridList_small)

# need this predict function for non-stationary, anisotropic awght
predict.multivariateSurfaceGrid <- function(object, x) {
  dimZ <- dim(object$z)
  L <- dimZ[3]
  out <- matrix(NA, nrow = nrow(x), ncol = L)
  for (l in 1:L) {
    out[, l] <- interp.surface(
      list(x = object$x, y = object$y, z = object$z[,, l]), x
    )
  }
  return(out)
}


# ------------------------------------
# Data & stationary LKinfo (outside loop)
# ------------------------------------
stations_df <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/eea_df_2023.rds")
cams_df     <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/cams_df_2023.rds")

# covariates used
covars <- c("rh", "ssr", "t2m",
            "windspeed", "sl_blh", "EM_NO2", "lag_cams_no2")

# results storage: test RMSE only
results <- data.frame(
  day                = 1:365,
  rmse_test_FRK   = NA_real_,
  time_sec           = NA_real_
)

# ------------------------------------
# Main loop over days
# ------------------------------------
for (i in 1:365) {
  t_start <- proc.time()
  chosen_day <- i
  
  cat("=== Day", chosen_day, "===\n")
  
  # SETUP FRK MODEL HERE
  
  
  
  # subset down to our chosen day 
  stations_df_day <- stations_df[
    stations_df$time == (stations_df$time[1] + chosen_day - 1),
  ]
  
  cams_df_day <- cams_df[
    cams_df$time == (cams_df$time[1] + chosen_day - 1),
  ]
  
  # if no CAMS data that day, skip
  if (nrow(cams_df_day) == 0) {
    cat("No CAMS data for day", chosen_day, "- skipping\n")
    results$time_sec[i] <- (proc.time() - t_start)[["elapsed"]]
    next
  }
  
  # scale ssr
  cams_df_day[, "ssr"]     <- cams_df_day[, "ssr"] / 10^6
  stations_df_day[, "ssr"] <- stations_df_day[, "ssr"] / 10^6
  
  # build sf objects
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
  train_idx   <- unique(nearest_idx)
  
  train_df <- cams_df_day[train_idx, ]
  test_df  <- cams_df_day[-train_idx, ]
  
  # if too few training points, skip day
  if (nrow(train_df) < 10) {
    cat("Too few training points for day", chosen_day, "- skipping\n")
    results$time_sec[i] <- (proc.time() - t_start)[["elapsed"]]
    next
  }
  
  # TRAIN
  s <- train_df[, c("Longitude", "Latitude")]
  y <- train_df$NO2
  Z <- as.matrix(train_df[, covars])
  
  
  
  # FIT FRK MODEL TO THE TRAINING DATA S Y Z
  
  
  # TEST
  s_test <- test_df[, c("Longitude", "Latitude")]
  y_test <- test_df$NO2
  Z_test <- as.matrix(test_df[, covars])
  
  
  
  # PREDICT FRK MODEL ON THE TEST SET 
  
  
  
  # CALCULATE FRK RMSE
  #rmse_test_FRK  <- sqrt(mean((FRK_PREDICTION_HERE  - y_test)^2))
  
  
  
  # timing
  elapsed <- (proc.time() - t_start)[["elapsed"]]
  
  # store results
  results[i, ] <- list(
    day                = chosen_day,
    rmse_test_FRK   = rmse_test_FRK,
    time_sec           = elapsed
  )
  
  cat("Finished day", chosen_day, "in", round(elapsed, 2), "seconds\n")
  print(results)
}



results_filtered <- results[!is.na(results$rmse_test_stat), ]
results <- results_filtered

# save if you want
saveRDS(results, "WHEREVERYOUWANTTOSAVETHEM")

