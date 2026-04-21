# Libraries
library(rhdf5)
library(LatticeKrig)
library(sf)
library(tictoc)

# for parallel option
library(foreach)
library(doParallel)

# ------------------------------------
# User choices
# ------------------------------------
n_days  <- 3    # how many sequential days to run
cv_type <- 10    # K-fold value
# ------------------------------------

set.seed(777)

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


# ---------------------------------------
# Function to estimate additional kappa2
# for an already nonstationary LKinfo
# ---------------------------------------

est_additional_kappa2 <- function(
    s, 
    y, 
    Z, 
    LKinfo, 
    kappa2_weights,
    verbose = FALSE
){

  kappa2_weights <- kappa2_weights
  
  LKinfo_additional <- LKinfo
  LKinfo_additional$a.wght <- NA
  LKinfo_additional$a.wghtObject <- NULL
  LKinfo_additional$a.wght <- list()
  LKinfo_additional$a.wght[[1]] <- 4.01
  attr(LKinfo_additional$a.wght, "isotropic") <- TRUE
  
  Lmodel_additional <- LatticeKrig_hacked(
    s,
    y,
    Z,
    LKinfo = LKinfo_additional, 
    findAwght = T, 
    verbose = verbose
  )
  
  LKinfo_new <- LKinfo
  LKinfo_new$a.wght[[1]][,5] <- LKinfo_new$a.wght[[1]][,5] + 
    (Lmodel_additional$LKinfo$a.wght[[1]]-4)
  
  return(LKinfo_new)
}


# ------------------------------------
# Load Initial Data 
# ------------------------------------
stations_df <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/eea_df_2023.rds")
stations_df$ssr <- stations_df$ssr / 10^6

covars <- c(
  "rh", "ssr", "t2m", 
  "windspeed", "sl_blh", "EM_NO2", 
  "DEM", "lag_cams_no2"
)

# for weights 
cams_df <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/cams_df_2023padded.rds")
cams_df <- cams_df[cams_df$time >= 19358 & cams_df$time <= 19722, ]
cams_df[,"ssr"] <- cams_df[,"ssr"] / 10^6
elev_map <-cams_df$DEM[cams_df$time == 19358]
elev_map_mask <- elev_map 
elev_map_mask[elev_map_mask != 0] <- 1
length(elev_map_mask)

par(mfrow = c(1,2))
image.plot(as.surface(sGrid_small, elev_map), col = magma(256))
image.plot(as.surface(sGrid_small, elev_map_mask), col = magma(256))
par(mfrow = c(1,1))

# ------------------------------------
# Experiment Loop 
# ------------------------------------

cv_results_df <- data.frame(
  day = 1:n_days,
  lm_rmse = NA_real_,
  stat_rmse = NA_real_,
  ns_rmse = NA_real_,
  ns_k2add_rmse = NA_real_
)

for (day_idx in 1:n_days) {
  cat("Starting day", day_idx, "of", n_days, "\n")
  tic()
  
  stations_df_day <- stations_df[
    stations_df$time == (stations_df$time[1] + day_idx - 1),
  ]
  # remove all rows where stations_df_day$EEA_NO2 is NA
  stations_df_day <- stations_df_day[!is.na(stations_df_day$EEA_NO2), ]
  
  n_obs <- nrow(stations_df_day)
  if (n_obs < cv_type * 2) {
    warning(paste("Day", day_idx, "has too few observations:", n_obs))
    next
  }
  
  # assign fold ids
  fold_ids <- sample(rep(1:cv_type, length.out = n_obs))
  
  # ------------------------------------
  # Model setup 
  # ------------------------------------
  
  # base, stationary comparison model
  LKinfo_stat <- LKrigSetup(
    sGrid_small, 
    a.wght = 4.01,
    NC = 128, 
    nlevel = 1,
    normalize = TRUE,
    NC.buffer = 0
  )
  
  # nonstationary model straight from STUN CTM param estimates
  LKinfo_ns <- make_nonstat_LKinfo(
    file_path   = "Italy_AQ_AmortisedLatticeKrig/data/STUN_param_df_2023padded.h5",
    dataset_name= "arx1_surround_30rep_output",
    day         = day_idx,
    gridlist    = gridList_small,
    normalize   = TRUE
  )
  nonstat_Q <- LKrig.precision(LKinfo_ns)
  
  # vectors to store RMSEs for this day across folds
  rmse_lm_vec       <- numeric(cv_type)
  rmse_stat_vec     <- numeric(cv_type)
  rmse_ns_vec       <- numeric(cv_type)
  rmse_ns_k2add_vec <- numeric(cv_type)
  
  # ------------------------------------
  # K-fold CV within this day
  # ------------------------------------
  for (k in 1:cv_type) {
    # setting up training and testing datasets and covars 
    test_ind <- which(fold_ids == k)
    train_df <- stations_df_day[-test_ind, ]
    test_df  <- stations_df_day[test_ind, ]
    
    s_train <- train_df[, c("Longitude", "Latitude")]
    y_train <- train_df$EEA_NO2
    Z_train <- as.matrix(train_df[, covars])
    
    s_test <- test_df[, c("Longitude", "Latitude")]
    y_test <- test_df$EEA_NO2
    Z_test <- as.matrix(test_df[, covars])
    
    # ------------------------------------
    # Fitting models 
    # ------------------------------------
    
    # linear model
    linear_model <- lm(
      EEA_NO2 ~ Longitude + Latitude + rh + ssr + t2m +
        windspeed + sl_blh + EM_NO2 + DEM + lag_cams_no2,
      data = train_df
    )
    
    # stationary LKrig with your custom LKinfo_stat
    Lmodel_stat <- LatticeKrig(
      s_train,
      y_train,
      Z = Z_train,
      LKinfo = LKinfo_stat,
      findAwght = TRUE
    )
    
    # original STUN nonstationary
    Lmodel_ns <- LatticeKrig(
      x = s_train,
      y = y_train,
      Z = Z_train,
      LKinfo = LKinfo_ns
    )
    
    # nonstat model with additonal kappa2 estimated 
    # default kappa2 estimation 
    # kappa2_weights <- 1
    
    # can also use the digital elevation mask 
    kappa2_weights <- elev_map_mask
    
    LKinfo_ns_k2add <- est_additional_kappa2(
      s_train,
      y_train,
      Z_train,
      LKinfo_ns, 
      kappa2_weights, 
      verbose = FALSE
    )
    print("Finished estimating additional kappa2")
    
    Lmodel_ns_k2add <- LatticeKrig(
      x = s_train,
      y = y_train,
      Z = Z_train,
      LKinfo = LKinfo_ns_k2add
    )
    
    # ------------------------------------
    # Evaluating models on test set
    # ------------------------------------
    
    pred_linear     <- predict(linear_model,     newdata = test_df)
    pred_stat       <- predict(Lmodel_stat,      s_test, Z = Z_test)
    pred_ns         <- predict(Lmodel_ns,        s_test, Z = Z_test)
    pred_ns_k2add   <- predict(Lmodel_ns_k2add,  s_test, Z = Z_test)
    
    rmse_lm_vec[k]       <- sqrt(mean((pred_linear   - y_test)^2))
    rmse_stat_vec[k]     <- sqrt(mean((pred_stat     - y_test)^2))
    rmse_ns_vec[k]       <- sqrt(mean((pred_ns       - y_test)^2))
    rmse_ns_k2add_vec[k] <- sqrt(mean((pred_ns_k2add - y_test)^2))
  }
  
  # ------------------------------------
  # Store mean RMSEs for this day
  # ------------------------------------
  cv_results_df$lm_rmse[day_idx]        <- mean(rmse_lm_vec,       na.rm = TRUE)
  cv_results_df$stat_rmse[day_idx]      <- mean(rmse_stat_vec,     na.rm = TRUE)
  cv_results_df$ns_rmse[day_idx]        <- mean(rmse_ns_vec,       na.rm = TRUE)
  cv_results_df$ns_k2add_rmse[day_idx]  <- mean(rmse_ns_k2add_vec, na.rm = TRUE)
  
  toc()
}


# ------------------------------------
# Save results
# ------------------------------------
saveRDS(cv_results_df, "Italy_AQ_AmortisedLatticeKrig/results/cv_results_df.rds")
print("Cross-validation finished. Results saved.")
print(cv_results_df)