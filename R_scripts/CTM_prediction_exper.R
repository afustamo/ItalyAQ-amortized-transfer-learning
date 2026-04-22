##########             IMPORTANT              ###########
#-------------------------------------------------------#
# Make sure to set your working directory 
# "To Source File Location". 
# Also, make sure you have downloaded all required 
# packages from required_packages.R. 
# Now you can run the script. 
#-------------------------------------------------------#

# Libraries
library(rhdf5)
library(LatticeKrig)
library(sf)
library(fields)
library(here)

# import helpful funcs
source(here("R_scripts", "helpful_functions.R"))

# seed for reproducibility 
set.seed(7)

# create grid 
rows_small <- length(seq(6.1, 18.8, by = .1))

gridList_small <- list(
  x = seq(6.1, 18.8, by = .1),
  y = seq(35.2, 47.9, by = .1)
)
sGrid_small <- make.surface.grid(gridList_small)


# ------------------------------------
# Data & stationary LKinfo (outside loop)
# ------------------------------------
stations_df <- readRDS(here("data", "eea_df_2023.rds"))
cams_df     <- readRDS(here("data", "cams_df_2023padded.rds"))

cams_df <- cams_df[cams_df$time >= 19358 & cams_df$time <= 19722, ]

# stationary LKinfo does not depend on day
LKinfo_stat <- LKrigSetup(
  sGrid_small,
  NC = 128,
  nlevel = 1,
  normalize = TRUE,
  NC.buffer = 0,
  a.wght = 4.1
)

# covariates used
covars <- c("rh", "ssr", "t2m",
            "windspeed", "sl_blh", "DEM",
            "EM_NO2", "lag_cams_no2")


# results storage: test RMSE only
results <- data.frame(
  day                = 1:365,
  rmse_test_linear   = NA_real_,
  rmse_test_LKdef    = NA_real_,
  rmse_test_ked      = NA_real_,
  rmse_test_stat     = NA_real_,
  rmse_test_ns       = NA_real_,   
  time_sec           = NA_real_, 
  time_coef_stat     = NA_real_,
  time_coef_ns       = NA_real_, 
  kappa2_stat         = NA_real_
)

# # results storage: full field predictions
predictions_df <- cams_df[,c("Longitude", "Latitude", "time","NO2")]
predictions_df$linear_pred <- NA_real_
predictions_df$LKdef_pred <- NA_real_
predictions_df$ked_pred <- NA_real_
predictions_df$stat_pred <- NA_real_
predictions_df$ns_pred <- NA_real_


# ------------------------------------
# Main loop over days
# ------------------------------------
for (i in c(1:365)) {
  t_start <- proc.time()
  chosen_day <- i
  
  cat("=== Day", chosen_day, "===\n")
  
  # ----- Nonstationary LKinfo: original STUN ns -----
  LKinfo_ns <- make_nonstat_LKinfo(
    file_path   = here("results", "STUN_param_df_2023padded.h5"),
    dataset_name= "arx1_surround_30rep_output",
    day         = chosen_day,
    gridlist    = gridList_small,
    normalize   = TRUE,
    sanity_plotting = FALSE,
    sanity_sim      = FALSE
  )

  
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
  
  # linear model
  linear_model <- lm(
    NO2 ~ Longitude + Latitude + rh + ssr + t2m +
      windspeed + sl_blh + EM_NO2 + DEM + lag_cams_no2,
    data = train_df
  )

  #LKdef 
  Lmodel_def <- LatticeKrig(
    s, 
    y, 
    Z = Z
    # findAwght = TRUE
  )

  # ked
  Lmodel_ked <- spatialProcess(
    s,
    y, 
    Z = Z
  )
  
  # stationary LKrig with your custom LKinfo_stat
  Lmodel_stat <- LatticeKrig(
    s,
    y,
    Z = Z,
    LKinfo = LKinfo_stat,
    findAwght = TRUE
  )
  
  # original STUN nonstationary
  Lmodel_ns <- LatticeKrig(
    x = s,
    y = y,
    Z = Z,
    LKinfo = LKinfo_ns
  )
  
  # TEST
  s_test <- test_df[, c("Longitude", "Latitude")]
  y_test <- test_df$NO2
  Z_test <- as.matrix(test_df[, covars])
  
  pred_linear     <- predict(linear_model,     newdata = test_df)
  pred_LKdef      <- predict(Lmodel_def,      s_test, Z = Z_test)
  pred_ked        <- predict(Lmodel_ked,      s_test, Z = Z_test)
  pred_stat       <- predict(Lmodel_stat,      s_test, Z = Z_test)
  pred_ns     <- predict(Lmodel_ns,    s_test, Z = Z_test)
  
  rmse_test_linear   <- sqrt(mean((pred_linear   - y_test)^2))
  rmse_test_LKdef    <- sqrt(mean((pred_LKdef    - y_test)^2))
  rmse_test_ked      <- sqrt(mean((pred_ked      - y_test)^2))
  rmse_test_stat     <- sqrt(mean((pred_stat     - y_test)^2))
  rmse_test_ns   <- sqrt(mean((pred_ns   - y_test)^2))

  kappa2_stat <- Lmodel_stat$LKinfo$a.wght[[1]]
  time_coef_stat <- Lmodel_stat$d.coef[11,1]
  time_coef_ns   <- Lmodel_ns$d.coef[11,1]
  
  # # timing
  # elapsed <- (proc.time() - t_start)[["elapsed"]]
  
  # store results
  results[i, ] <- list(
    day                = chosen_day,
    rmse_test_linear   = rmse_test_linear,
    rmse_test_LKdef    = rmse_test_LKdef,
    rmse_test_ked      = rmse_test_ked,
    rmse_test_stat     = rmse_test_stat,
    rmse_test_ns       = rmse_test_ns,
    time_sec           = elapsed, 
    time_coef_stat     = time_coef_stat,
    time_coef_ns       = time_coef_ns,
    kappa2_stat        = kappa2_stat
  )
  
  
  s_full <- cams_df_day[, c("Longitude", "Latitude")]
  Z_full <- as.matrix(cams_df_day[, covars])
  
  # # full preds
  pred_linear_full  <- predict(linear_model,     newdata = cams_df_day)
  pred_LKdef_full   <- predict(Lmodel_def,      s_full, Z = Z_full)
  pred_ked_full     <- predict(Lmodel_ked,      s_full, Z = Z_full)
  pred_stat_full    <- predict(Lmodel_stat,      s_full, Z = Z_full)
  pred_ns_full      <- predict(Lmodel_ns,        s_full, Z = Z_full)
  
  
  # # time adjust 
  current_time <- cams_df_day$time[1]
  day_mask <- predictions_df$time == current_time
  
  # # save full preds 
  predictions_df$linear_pred[day_mask]  <- pred_linear_full
  predictions_df$LKdef_pred[day_mask]   <- pred_LKdef_full
  predictions_df$ked_pred[day_mask]     <- pred_ked_full
  predictions_df$stat_pred[day_mask] <- pred_stat_full
  predictions_df$ns_pred[day_mask]   <- pred_ns_full
  

  # timing
  elapsed <- (proc.time() - t_start)[["elapsed"]]
  
  cat("Finished day", chosen_day, "in", round(elapsed, 2), "seconds\n")
  # print(results)
}

# save if you want
saveRDS(results, here("results", "rmse_CTM_reconstruct.rds"))
saveRDS(predictions_df, here("results", "predictions_CTM_reconstruct.rds"))