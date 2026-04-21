# Libraries
library(rhdf5)
library(LatticeKrig)
library(sf)
library(fields)

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
# Function to create an LKinfo object
# from STUN / CNN parameters
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
  params <- params[, 128:1, , day]
  
  # recover kappa
  kappa2 <- exp(params[,,1])
  # just in case we need awght
  awght <- kappa2 + 4
  # theta needs to be transformed like this 
  theta <- params[,,2] + pi/2
  rho   <- params[,,3]
  
  if (sanity_plotting){
    par(mfrow = c(1,3))
    imagePlot(as.surface(gridlist, kappa2),          main = "kappa2",     col = viridis(256))
    world(add = TRUE, col = "white", lwd = 1)
    imagePlot(as.surface(gridlist, theta - pi/2),    main = "theta(adj)", col = viridis(256))
    world(add = TRUE, col = "white", lwd = 1)
    imagePlot(as.surface(gridlist, rho),             main = "rho",        col = viridis(256))
    world(add = TRUE, col = "white", lwd = 1)
    par(mfrow = c(1,1))
  }
  
  # need these for encoding into LK
  rhox <- sqrt(rho)
  rhoy <- 1 / rhox
  
  # create H tensor out of params
  H11 <- (rhox^2 * (cos(theta))^2) + (rhoy^2 * (sin(theta))^2)
  H12 <- (rhoy^2 - rhox^2) * (sin(theta) * cos(theta))
  H21 <- H12 
  H22 <- (rhox^2 * (sin(theta))^2) + (rhoy^2 * (cos(theta))^2)
  
  rows <- length(gridlist$x)
  
  # fill the high dimensional stencil (9 fields)
  stencil_tensor <- array(NA, c(rows, rows, 9))
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
  awght_obj <- list(x = gridlist$x, y = gridlist$y, z = stencil_tensor)
  class(awght_obj) <- "multivariateSurfaceGrid"
  
  sGrid <- make.surface.grid(gridlist)
  
  LKinfo <- LKrigSetup(
    # dont change grid
    sGrid, 
    # you change awght indirectly with day and file selection
    a.wghtObject = awght_obj,
    # the rest of these you change directly in function calls 
    NC = NC, 
    nlevel = nlevel,
    normalize = normalize,
    NC.buffer = NC.buffer
  )
  
  if (sanity_sim){
    test <- LKrig.sim(
      sGrid,
      LKinfo = LKinfo,
      M = 1
    )
    
    if (sanity_plotting){
      imagePlot(as.surface(gridlist, test), col = turbo(256))
      world(add = TRUE, col = "black", lwd = 1)
    }
  }
  
  return(LKinfo)
}

# ------------------------------------
# Data & stationary LKinfo (outside loop)
# ------------------------------------
stations_df <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/eea_df_2023.rds")
cams_df     <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/cams_df_2023padded.rds")

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
  # rmse_test_ns_s1    = NA_real_,
  # rmse_test_ns_r1    = NA_real_,
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

# #new
# predictions_df$lm_SE <- NA_real_
# predictions_df$stat_SE <- NA_real_
# predictions_df$ns_SE <- NA_real_
# predictions_df$stat_upr <- NA_real_
# predictions_df$stat_lwr <- NA_real_
# predictions_df$ns_upr <- NA_real_
# predictions_df$ns_lwr <- NA_real_


# ------------------------------------
# Main loop over days
# ------------------------------------
for (i in c(1:365)) {
  t_start <- proc.time()
  chosen_day <- i
  
  cat("=== Day", chosen_day, "===\n")
  
  # ----- Nonstationary LKinfo: original STUN ns -----
  LKinfo_ns <- make_nonstat_LKinfo(
    file_path   = "Italy_AQ_AmortisedLatticeKrig/data/STUN_param_df_2023padded.h5",
    dataset_name= "arx1_surround_30rep_output",
    day         = chosen_day,
    gridlist    = gridList_small,
    normalize   = TRUE,
    sanity_plotting = FALSE,
    sanity_sim      = FALSE
  )

  # LKinfo_ns_s1 <- make_nonstat_LKinfo(
  #   file_path   = "Italy_AQ_AmortisedLatticeKrig/data/old/STUN_param_df_2023.h5",
  #   dataset_name= "arx1_surround_1rep_output",
  #   day         = chosen_day,
  #   gridlist    = gridList_small,
  #   normalize   = TRUE,
  #   sanity_plotting = FALSE,
  #   sanity_sim      = FALSE
  # )

  # LKinfo_ns_r1 <- make_nonstat_LKinfo(
  #   file_path   = "Italy_AQ_AmortisedLatticeKrig/data/old/STUN_param_df_2023.h5",
  #   dataset_name= "arx1_rscale_1rep_output",
  #   day         = chosen_day,
  #   gridlist    = gridList_small,
  #   normalize   = TRUE,
  #   sanity_plotting = FALSE,
  #   sanity_sim      = FALSE
  # )
  
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
  
  # timing
  elapsed <- (proc.time() - t_start)[["elapsed"]]
  
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
  
  # # standard error (NEW)
  # SE_linear_full <- predict(linear_model, newdata = cams_df_day, se.fit = TRUE)$se.fit
  # SE_stat_full   <- predictSE.LKrig(Lmodel_stat, s_full, Z = Z_full)
  # SE_ns_full     <- predictSE.LKrig(Lmodel_ns,   s_full, Z = Z_full)
  
  # stat_upper <- pred_stat_full + 1.96 * sqrt(SE_stat_full^2 + Lmodel_stat$tau.MLE^2)
  # stat_lower <- pred_stat_full - 1.96 * sqrt(SE_stat_full^2 + Lmodel_stat$tau.MLE^2)
  # ns_upper <- pred_ns_full + 1.96 * sqrt(SE_ns_full^2 + Lmodel_ns$tau.MLE^2)
  # ns_lower <- pred_ns_full - 1.96 * sqrt(SE_ns_full^2 + Lmodel_ns$tau.MLE^2)
  
  
  # # time adjust 
  current_time <- cams_df_day$time[1]
  day_mask <- predictions_df$time == current_time
  
  # # save full preds 
  predictions_df$linear_pred[day_mask]  <- pred_linear_full
  predictions_df$LKdef_pred[day_mask]   <- pred_LKdef_full
  predictions_df$ked_pred[day_mask]     <- pred_ked_full
  predictions_df$stat_pred[day_mask] <- pred_stat_full
  predictions_df$ns_pred[day_mask]   <- pred_ns_full
  
  # #NEW 
  # predictions_df$lm_SE[day_mask]       <- SE_linear_full
  # predictions_df$LK_stat_SE[day_mask]  <- SE_stat_full
  # predictions_df$LK_ns_SE[day_mask]    <- SE_ns_full
  # predictions_df$LK_stat_upr[day_mask] <- stat_upper
  # predictions_df$LK_stat_lwr[day_mask] <- stat_lower
  # predictions_df$LK_ns_upr[day_mask]   <- ns_upper
  # predictions_df$LK_ns_lwr[day_mask]   <- ns_lower
  
  # # # timing
  # elapsed <- (proc.time() - t_start)[["elapsed"]]
  
  cat("Finished day", chosen_day, "in", round(elapsed, 2), "seconds\n")
  # print(results)
}


# save if you want
saveRDS(results, "Italy_AQ_AmortisedLatticeKrig/results/rmse_CTM_reconstruct.rds")
saveRDS(predictions_df, "Italy_AQ_AmortisedLatticeKrig/results/predictions_CTM_reconstruct.rds")