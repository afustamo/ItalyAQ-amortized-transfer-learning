# Libraries
library(rhdf5)
library(LatticeKrig)
library(sf)

# import hacked LK files for additional kappa2 estimation
setwd("C:/Users/anton/Desktop/Research/Italy_AQ_paper")
source("Italy_AQ_AmortisedLatticeKrig/R_scripts/hacked_LatticeKrig.R")

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
# from STUN parameters
# ------------------------------------

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

# ------------------------------------
# Data & stationary LKinfo (outside loop)
# ------------------------------------
stations_df <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/eea_df_2023.rds")
cams_df     <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/cams_df_2023.rds")

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
            "windspeed", "sl_blh", "EM_NO2", "lag_cams_no2")

# results storage (now includes ns, ns1, ns2)
results <- data.frame(
  day = 1:365,
  rmse_test_linear = NA_real_,
  rmse_test_stat   = NA_real_,
  rmse_test_ns     = NA_real_,
  rmse_test_ns1    = NA_real_,
  rmse_test_ns2    = NA_real_,
  rmse_full_linear = NA_real_,
  rmse_full_stat   = NA_real_,
  rmse_full_ns     = NA_real_,
  rmse_full_ns1    = NA_real_,
  rmse_full_ns2    = NA_real_,
  time_sec         = NA_real_
)

# ------------------------------------
# Main loop over days
# ------------------------------------
for (i in sample(1:365, 100)) {
  t_start <- proc.time()
  chosen_day <- i
  
  cat("=== Day", chosen_day, "===\n")
  
  # Nonstationary LKinfo for this day: original model (ns)
  LKinfo_ns <- make_nonstat_LKinfo(
    file_path = "Italy_AQ_AmortisedLatticeKrig/data/STUN_param_df_2023.h5",
    dataset_name = "arx1_surround_30rep_output",
    day = chosen_day,
    gridlist = gridList_small,
    normalize = TRUE,
    sanity_plotting = FALSE,
    sanity_sim = FALSE
  )
  
  # Nonstationary LKinfo for this day: model 1
  LKinfo_ns1 <- make_nonstat_LKinfo(
    file_path = "Italy_AQ_AmortisedLatticeKrig/data/STUN_param1_df_2023.h5",
    dataset_name = "arx1_surround1_30rep_output",
    day = chosen_day,
    gridlist = gridList_small,
    normalize = TRUE,
    sanity_plotting = FALSE,
    sanity_sim = FALSE
  )
  
  # Nonstationary LKinfo for this day: model 2
  LKinfo_ns2 <- make_nonstat_LKinfo(
    file_path = "Italy_AQ_AmortisedLatticeKrig/data/STUN_param2_df_2023.h5",
    dataset_name = "arx1_surround2_30rep_output",
    day = chosen_day,
    gridlist = gridList_small,
    normalize = TRUE,
    sanity_plotting = FALSE,
    sanity_sim = FALSE
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
      windspeed + sl_blh + EM_NO2 + lag_cams_no2,
    data = train_df
  )
  
  # LK models
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
  
  Lmodel_ns1 <- LatticeKrig(
    x = s,
    y = y,
    Z = Z,
    LKinfo = LKinfo_ns1
  )
  
  Lmodel_ns2 <- LatticeKrig(
    x = s,
    y = y,
    Z = Z,
    LKinfo = LKinfo_ns2
  )
  
  # TEST
  s_test <- test_df[, c("Longitude", "Latitude")]
  y_test <- test_df$NO2
  Z_test <- as.matrix(test_df[, covars])
  
  pred_linear <- predict(linear_model, newdata = test_df)
  pred_stat   <- predict(Lmodel_stat, s_test, Z = Z_test)
  pred_ns     <- predict(Lmodel_ns,   s_test, Z = Z_test)
  pred_ns1    <- predict(Lmodel_ns1,  s_test, Z = Z_test)
  pred_ns2    <- predict(Lmodel_ns2,  s_test, Z = Z_test)
  
  rmse_test_linear <- sqrt(mean((pred_linear - y_test)^2))
  rmse_test_stat   <- sqrt(mean((pred_stat   - y_test)^2))
  rmse_test_ns     <- sqrt(mean((pred_ns     - y_test)^2))
  rmse_test_ns1    <- sqrt(mean((pred_ns1    - y_test)^2))
  rmse_test_ns2    <- sqrt(mean((pred_ns2    - y_test)^2))
  
  # FULL GRID
  Z_map     <- as.matrix(cams_df_day[, covars])
  cams_grid <- make.surface.grid(gridList_small)
  
  surf_pred_stat    <- predict(Lmodel_stat, cams_grid, Z = Z_map)
  surf_pred_ns      <- predict(Lmodel_ns,   cams_grid, Z = Z_map)
  surf_pred_ns1     <- predict(Lmodel_ns1,  cams_grid, Z = Z_map)
  surf_pred_ns2     <- predict(Lmodel_ns2,  cams_grid, Z = Z_map)
  surf_pred_linear  <- predict(linear_model, newdata = cams_df_day)
  
  rmse_full_stat    <- sqrt(mean((surf_pred_stat   - cams_df_day$NO2)^2))
  rmse_full_ns      <- sqrt(mean((surf_pred_ns     - cams_df_day$NO2)^2))
  rmse_full_ns1     <- sqrt(mean((surf_pred_ns1    - cams_df_day$NO2)^2))
  rmse_full_ns2     <- sqrt(mean((surf_pred_ns2    - cams_df_day$NO2)^2))
  rmse_full_linear  <- sqrt(mean((surf_pred_linear - cams_df_day$NO2)^2))
  
  # timing
  elapsed <- (proc.time() - t_start)[["elapsed"]]
  
  # store results
  results[i, ] <- list(
    day = chosen_day,
    rmse_test_linear = rmse_test_linear,
    rmse_test_stat   = rmse_test_stat,
    rmse_test_ns     = rmse_test_ns,
    rmse_test_ns1    = rmse_test_ns1,
    rmse_test_ns2    = rmse_test_ns2,
    rmse_full_linear = rmse_full_linear,
    rmse_full_stat   = rmse_full_stat,
    rmse_full_ns     = rmse_full_ns,
    rmse_full_ns1    = rmse_full_ns1,
    rmse_full_ns2    = rmse_full_ns2,
    time_sec         = elapsed
  )
  
  cat("Finished day", chosen_day, "in", round(elapsed, 2), "seconds\n")
  # print where results is not na  
  print(results[!is.na(results$rmse_test_stat), ])
}

results_filtered <- results[!is.na(results$rmse_test_stat), ]
results <- results_filtered

# save results if you want
# saveRDS(results, "Italy_AQ_AmortisedLatticeKrig/results/rmse_STUN_LK_365days_threeNS.rds")
# write.csv(results, "Italy_AQ_AmortisedLatticeKrig/results/rmse_STUN_LK_365days_threeNS.csv", row.names = FALSE)

## Example comparisons

# ns vs stat
nrow(results[results$rmse_full_ns < results$rmse_full_stat, ])
mean(results$rmse_full_ns,  na.rm = TRUE) / mean(results$rmse_full_stat, na.rm = TRUE)
t.test(results$rmse_full_ns, results$rmse_full_stat, paired = TRUE)

# ns1 vs ns
nrow(results[results$rmse_full_ns1 < results$rmse_full_ns, ])
mean(results$rmse_full_ns1, na.rm = TRUE) / mean(results$rmse_full_ns, na.rm = TRUE)
t.test(results$rmse_full_ns1, results$rmse_full_ns, paired = TRUE)

# ns2 vs ns1
nrow(results[results$rmse_full_ns2 < results$rmse_full_ns1, ])
mean(results$rmse_full_ns2, na.rm = TRUE) / mean(results$rmse_full_ns1, na.rm = TRUE)
t.test(results$rmse_full_ns2, results$rmse_full_ns1, paired = TRUE)
