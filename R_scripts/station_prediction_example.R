# Libraries
library(rhdf5)
library(LatticeKrig)
library(sf)
library(tictoc)
library(scico)

# ------------------------------------
# Choose which day
# ------------------------------------
chosen_day <- 305 #14 #305 is nov 1. 103 is apr 13

# import hacked LK files for additional kappa2 estimation
setwd("C:/Users/anton/Desktop/Research/Italy_AQ_paper")
source("Italy_AQ_AmortisedLatticeKrig/R_scripts/hacked_LatticeKrig.R")

set.seed(66)

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
  # assign("kappa2_weights", as.numeric(kappa2_weights), envir = .GlobalEnv)
  
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
    findAwght = T
  )
  
  LKinfo_new <- LKinfo
  # LKinfo_new$a.wght[[1]][,5] <- LKinfo_new$a.wght[[1]][,5] +
  #   (kappa2_weights * (Lmodel_additional$LKinfo$a.wght[[1]]-4))
  
  LKinfo_new$a.wght[[1]][,5] <- LKinfo_new$a.wght[[1]][,5] +
    (Lmodel_additional$LKinfo$a.wght[[1]]-4)
  
  return(LKinfo_new)
}


# ------------------------------------
# Dataset setup 
# ------------------------------------
stations_df <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/eea_df_2023.rds")
stations_df[, "ssr"] <- stations_df[, "ssr"] / 10^6

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

# subset to chosen day 
stations_df_day <- stations_df[
  stations_df$time == (stations_df$time[1] + chosen_day - 1),
]
# remove all rows where stations_df_day$EEA_NO2 is NA
stations_df_day <- stations_df_day[!is.na(stations_df_day$EEA_NO2), ]


# setting up training and testing datasets and covars 
covars <- c("rh", "ssr", "t2m",
            "windspeed", "sl_blh", 
            "DEM", "EM_NO2", "lag_cams_no2")

test_ind <- sample(1:nrow(stations_df_day), floor(nrow(stations_df_day)/10) )
train_df <- stations_df_day[-test_ind, ]
test_df  <- stations_df_day[test_ind, ]

s_train <- train_df[, c("Longitude", "Latitude")]
y_train <- train_df$EEA_NO2
Z_train <- as.matrix(train_df[, covars])

s_test <- test_df[, c("Longitude", "Latitude")]
y_test <- test_df$EEA_NO2
Z_test <- as.matrix(test_df[, covars])


# ------------------------------------
# Model setup 
# ------------------------------------


# base, stationary comparison model
LKinfo_stat <- LKrigSetup(
  sGrid_small, 
  a.wght = 4.01,
  # the rest of these you change directly in function calls 
  NC = 128, 
  nlevel = 1,
  normalize = TRUE,
  NC.buffer = 0
)


# nonstationary model straight from STUN CTM param estimates
LKinfo_ns <- make_nonstat_LKinfo(
  file_path   = "Italy_AQ_AmortisedLatticeKrig/data/STUN_param_df_2023padded.h5",
  dataset_name= "arx1_surround_30rep_output",
  day         = chosen_day,
  gridlist    = gridList_small,
  normalize   = TRUE,
  sanity_plotting = FALSE,
  sanity_sim      = FALSE
)
nonstat_Q <- LKrig.precision(LKinfo_ns)


# nonstat model with additonal kappa2 estimated 
#default
# kappa2_weights <- 1

# land ocean
kappa2_weights <- elev_map_mask

LKinfo_ns_k2add <- est_additional_kappa2(
  s_train,
  y_train,
  Z_train,
  LKinfo_ns, 
  kappa2_weights = kappa2_weights
)

# checking additional kappa2 for sanity
for (i in 1:9){
  test.for.zero(LKinfo_ns_k2add$a.wght[[1]][,i], LKinfo_ns$a.wght[[1]][,i])
}

print(LKinfo_ns_k2add$a.wght[[1]][,5] - LKinfo_ns$a.wght[[1]][,5])
min(LKinfo_ns_k2add$a.wght[[1]][,5] - LKinfo_ns$a.wght[[1]][,5])
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

Lmodel_ns_k2add <- LatticeKrig(
  x = s_train,
  y = y_train,
  Z = Z_train,
  LKinfo = LKinfo_ns_k2add
)

# summary of fits on train data 
summary(abs(linear_model$residuals))
summary(abs(Lmodel_stat$residuals))
summary(abs(Lmodel_ns$residuals))
summary(abs(Lmodel_ns_k2add$residuals))



test1 <- Lmodel_ns$LKinfo$a.wghtObject$z[,,5] + 
  (2 * Lmodel_ns$LKinfo$a.wghtObject$z[,,2]) + 
  (2 * Lmodel_ns$LKinfo$a.wghtObject$z[,,4])

test2 <- Lmodel_ns_k2add$LKinfo$a.wghtObject$z[,,5] + 
  (2 * Lmodel_ns_k2add$LKinfo$a.wghtObject$z[,,2]) + 
  (2 * Lmodel_ns_k2add$LKinfo$a.wghtObject$z[,,4])

test11 <- Lmodel_ns$LKinfo$a.wght[[1]][,5] + 
  (2 * Lmodel_ns$LKinfo$a.wght[[1]][,2]) + 
  (2 * Lmodel_ns$LKinfo$a.wght[[1]][,4])

test22 <- Lmodel_ns_k2add$LKinfo$a.wght[[1]][,5] + 
  (2 * Lmodel_ns_k2add$LKinfo$a.wght[[1]][,2]) + 
  (2 * Lmodel_ns_k2add$LKinfo$a.wght[[1]][,4])

image.plot(as.surface(gridList_small, test11), col = viridis(256), main = "ns")
image.plot(as.surface(gridList_small, test22), col = viridis(256), main = "ns_k2add")

test3 <- Lmodel_stat$LKinfo$a.wght[[1]]-4

summary(test22 - test3)

image.plot(as.surface(gridList_small, test22 -test3), col = viridis(256), main = "ns_k2add")

# ------------------------------------
# Evaluating models on test set
# ------------------------------------

pred_linear     <- predict(linear_model,     newdata = test_df)
pred_stat       <- predict(Lmodel_stat,      s_test, Z = Z_test)
pred_ns         <- predict(Lmodel_ns,        s_test, Z = Z_test)
pred_ns_k2add  <- predict(Lmodel_ns_k2add,   s_test, Z = Z_test)

rmse_test_linear   <- sqrt(mean((pred_linear   - y_test)^2))
rmse_test_stat     <- sqrt(mean((pred_stat     - y_test)^2))
rmse_test_ns       <- sqrt(mean((pred_ns       - y_test)^2))
rmse_test_ns_k2add <- sqrt(mean((pred_ns_k2add - y_test)^2))

print(
  paste0(
    "Lm: ", round(rmse_test_linear, 4), "  ",
    "Stat: ", round(rmse_test_stat, 4), "  ",
    "Nonstat: ", round(rmse_test_ns, 4), "  ",
    "Nonstat additional kappa2: ", round(rmse_test_ns_k2add, 4)
  )
)


# ----------------------------------------
# Downscaled, gridded prediction surfaces
# ----------------------------------------

# loading in data
preds_df <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/pred_df_2023.rds")
preds_df[,"ssr"] <- preds_df[,"ssr"] / 10^6

# setup
preds_df_day <- preds_df[
  preds_df$time == (preds_df$time[1] + chosen_day - 1),
]

s_pred <- preds_df_day[, c("Longitude", "Latitude")]
Z_pred <- as.matrix(preds_df_day[, covars])

gridList_big <- list( 
  x = seq( 6.1,18.8,by=.05), 
  y= seq( 35.2,47.9,by=.05)
)

sGrid_big<- make.surface.grid(gridList_big)

surf_lm <- predict(
  linear_model, 
  newdata = preds_df_day
)

# surface predictions
tic()
surf_stat <- predict(Lmodel_stat, s_pred, Z = Z_pred)
toc()
surf_stat_fxd <- predict(Lmodel_stat, s_pred, Z_pred, just.fixed = TRUE)
surf_stat_gp <- surf_stat - surf_stat_fxd

tic()
surf_ns <- predict(Lmodel_ns, s_pred, Z = Z_pred)
toc()
surf_ns_fxd <- predict(Lmodel_ns, s_pred, Z_pred, just.fixed = TRUE)
surf_ns_gp <- surf_ns - surf_ns_fxd

tic()
surf_ns_k2add <- predict(Lmodel_ns_k2add, s_pred, Z = Z_pred)
toc()
surf_ns_k2add_fxd <- predict(Lmodel_ns_k2add, s_pred, Z_pred, just.fixed = TRUE)
surf_ns_k2add_gp <- surf_ns_k2add - surf_ns_k2add_fxd

# set zlimits for full surfaces
pred_max <- max(
  c(surf_lm, surf_stat, surf_ns, surf_ns_k2add),
  na.rm = TRUE
)
pred_min <- min(
  c(surf_lm, surf_stat, surf_ns, surf_ns_k2add),
  na.rm = TRUE
)

pred_fxd_max <- max(
  c(surf_lm, surf_stat_fxd, surf_ns_fxd, surf_ns_k2add_fxd),
  na.rm = TRUE
)
pred_fxd_min <- min(
  c(surf_lm, surf_stat_fxd, surf_ns_fxd, surf_ns_k2add_fxd),
  na.rm = TRUE
)

pred_gp_max <- max(
  c(surf_stat_gp, surf_ns_gp, surf_ns_k2add_gp),
  na.rm = TRUE
)
pred_gp_min <- min(
  c(surf_stat_gp, surf_ns_gp, surf_ns_k2add_gp),
  na.rm = TRUE
)


# full prediction surfaces
par(mfrow = c(1,4))
par(mar = c(2,2,2,2))
imagePlot(
  as.surface(gridList_big, surf_lm), 
  col=plasma(256), main = "LM Pred", 
  zlim = c(pred_min, pred_max)
)
world(add = TRUE, col = "black", lwd = 1)
imagePlot(
  as.surface(gridList_big, surf_stat),
  col=plasma(256), main = "Stat Pred",
  zlim = c(pred_min, pred_max)
)
world(add = TRUE, col = "black", lwd = 1)
imagePlot(
  as.surface(gridList_big, surf_ns),
  col = plasma(256), main = "Nonstat Pred",
  zlim = c(pred_min, pred_max)
)
world(add = TRUE, col = "black", lwd = 1)
imagePlot(
  as.surface(gridList_big, surf_ns_k2add),
  col = plasma(256), main = "Nonstat + kappa2 Pred",
  zlim = c(pred_min, pred_max)
)
world(add = TRUE, col = "black", lwd = 1)
par(mfrow = c(1,1))




# full preds, fixed effects, and gps only mega-image
par(mfrow = c(3,3))
# full preds 
imagePlot(
  as.surface(gridList_big, surf_stat),
  col=plasma(256), main = "Stat Pred",
  zlim = c(pred_min, pred_max)
)
world(add = TRUE, col = "black", lwd = 1)
imagePlot(
  as.surface(gridList_big, surf_ns),
  col = plasma(256), main = "Nonstat Pred",
  zlim = c(pred_min, pred_max)
)
world(add = TRUE, col = "black", lwd = 1)
imagePlot(
  as.surface(gridList_big, surf_ns_k2add),
  col = plasma(256), main = "Nonstat + kappa2 Pred",
  zlim = c(pred_min, pred_max)
)
world(add = TRUE, col = "black", lwd = 1)

# fixed components
imagePlot(
  as.surface(gridList_big, surf_lm), 
  col=inferno(256), main = "Stat Fixed",
  zlim = c(pred_min, pred_max)
)
world(add = TRUE, col = "black", lwd = 1)

imagePlot(
  as.surface(gridList_big, surf_ns_fxd),
  col=inferno(256), main = "Nonstat Fixed",
  zlim = c(pred_min, pred_max)
)
world(add = TRUE, col = "black", lwd = 1)

imagePlot(
  as.surface(gridList_big, surf_ns_k2add_fxd),
  col=inferno(256), main = "Nonstat + kappa2 Fixed",
  zlim = c(pred_min, pred_max)
)
world(add = TRUE, col = "black", lwd = 1)

# gp only 
imagePlot(
  as.surface(gridList_big, surf_stat_gp),
  col = rocket(256), main = "Stat GP",
  zlim = c(pred_gp_min, pred_gp_max)
)
world(add = TRUE, col = "black", lwd = 1)
imagePlot(
  as.surface(gridList_big, surf_ns_gp),
  col = rocket(256), main = "Nonstat GP",
  zlim = c(pred_gp_min, pred_gp_max)
)
world(add = TRUE, col = "black", lwd = 1)
imagePlot(
  as.surface(gridList_big, surf_ns_k2add_gp),
  col = rocket(256), main = "Nonstat + kappa2 GP",
  zlim = c(pred_gp_min, pred_gp_max)
)
world(add = TRUE, col = "black", lwd = 1)
par(mfrow = c(1,1))



# just some differences between final ns and stat
par(mfrow = c(1,2))
image.plot(
  as.surface(gridList_big, surf_ns_k2add_gp - surf_stat_gp),
  col = scico(256, palette = 'vik'), main = "Ns_k2add - Stat (GP only)"
)

image.plot(
  as.surface(gridList_big, surf_ns_k2add - surf_stat),
  col = scico(256, palette = 'vik'), main = "Ns_k2add - Stat (full preds)"
)
par(mfrow = c(1,1))


# ----------------------------------------
# Conditional Sim for SE surfaces
# ----------------------------------------

n_csims <- 1000 #1000

csim_stat <- LKrig.sim.conditional(
  Lmodel_stat, 
  M = n_csims, 
  x.grid = sGrid_big, 
  Z.grid = Z_pred
)
csim_ns <- LKrig.sim.conditional(
  Lmodel_ns, 
  M = n_csims, 
  x.grid = sGrid_big, 
  Z.grid = Z_pred
)
tic()
csim_ns_k2add <- LKrig.sim.conditional(
  Lmodel_ns_k2add, 
  M = n_csims, 
  x.grid = sGrid_big, 
  Z.grid = Z_pred
)
toc()

chosen_draw <- 7

par(mfrow = c(3,3))
image.plot(
  as.surface(sGrid_big, csim_stat$ghat), 
  col = plasma(256), main = "Stat ghat"
)
image.plot(
  as.surface(sGrid_big, csim_stat$g.draw[,chosen_draw]), 
  col = plasma(256), main = "Stat Sample Draw"
)
image.plot(
  as.surface(sGrid_big, csim_stat$SE), 
  col = mako(256), main = "Stat SE"
)

image.plot(
  as.surface(sGrid_big, csim_ns$ghat), 
  col = plasma(256), main = "Nonstat ghat"
)
image.plot(
  as.surface(sGrid_big, csim_ns$g.draw[,chosen_draw]), 
  col = plasma(256), main = "Nonstat Sample Draw"
)
image.plot(
  as.surface(sGrid_big, csim_ns$SE), 
  col = mako(256), main = "Nonstat SE"
)

image.plot(
  as.surface(sGrid_big, csim_ns_k2add$ghat), 
  col = plasma(256), main = "Nonstat + k2add ghat"
)
image.plot(
  as.surface(sGrid_big, csim_ns_k2add$g.draw[,chosen_draw]), 
  col = plasma(256), main = "Nonstat + k2add Sample Draw"
)
image.plot(
  as.surface(sGrid_big, csim_ns_k2add$SE), 
  col = mako(256), main = "Nonstat + k2add SE"
)
par(mfrow = c(1,1))





library(cmocean)
library(maps)

par(mfrow = c(1,4))

zlim_test <- c(
  min(surf_ns_k2add, stations_df_day$EEA_NO2),
  max(surf_ns_k2add, stations_df_day$EEA_NO2)
)

imagePlot(
  as.surface(gridList_big, surf_ns_k2add),
  col = plasma(256), main = "Full Pred (also ghat csim)", 
  zlim = zlim_test
)
map("world", add = TRUE, lwd = 1.2, col = "white")

# points(
#   x = stations_df_day$Longitude,
#   y = stations_df_day$Latitude,
#   pch = 1,        # "x"
#   col = "white",
#   cex = 0.4,
#   lwd = 0.5
# )

points(
  x = stations_df_day$Longitude,
  y = stations_df_day$Latitude,
  bg = stations_df_day$EEA_NO2,
  pch = 19,        # "x"
  col = plasma(256),
  cex = 0.6,
  lwd = 0.5
)

imagePlot(
  as.surface(gridList_big, surf_ns_k2add),
  col = plasma(256), main = "Full Pred (also ghat csim)"
)
map("world", add = TRUE, lwd = 1.2, col = "white")

points(
  x = stations_df_day$Longitude,
  y = stations_df_day$Latitude,
  pch = 1,        # "x"
  col = "white",
  cex = 0.4,
  lwd = 0.5
)


# image.plot(
#   as.surface(sGrid_big, csim_ns_k2add$ghat), 
#   col = plasma(256), main = "Ghat (csim)"
# )
# map("world", add = TRUE, lwd = 1.2, col = "white")

image.plot(
  as.surface(sGrid_big, csim_ns_k2add$g.draw[,chosen_draw]), 
  col = plasma(256), main = "Sample Draw (csim)"
)
map("world", add = TRUE, lwd = 1.2, col = "white")

points(
  x = stations_df_day$Longitude,
  y = stations_df_day$Latitude,
  pch = 1,        # "x"
  col = "white",
  cex = 0.4,
  lwd = 0.5
)

se_col <- rev(cmocean('deep')(256))
# se_col <- cmocean('haline')(256)
# se_col <- mako(256)

image.plot(
  as.surface(sGrid_big, csim_ns_k2add$SE), 
  col = se_col, main = "SE (csim)"
)
map("world", add = TRUE, lwd = 1.2, col = "white")

# white dots at station locations
points(
  x = stations_df_day$Longitude,
  y = stations_df_day$Latitude,
  pch = 1,        # "x"
  col = "white",
  cex = 0.4,
  lwd = 0.5
)

par(mfrow = c(1,1))








summary(csim_stat$SE)
summary(csim_ns_k2add$SE)
summary(csim_ns$SE)

image.plot(
  as.surface(sGrid_big, csim_ns_k2add$SE - csim_ns$SE), 
  col = se_col, main = "SE (csim)"
)
map("world", add = TRUE, lwd = 1.2, col = "white")


pdf("Paper/fine_grid_pred_nov1.pdf", width = 7, height = 5)
par(
  mfrow = c(1,2), 
  mar = c(4,3.6,2,2.2), 
  oma = c(0,0,2.0,0), 
  mgp   = c(2, 0.6, 0),
  bg = "white"
)

imagePlot(
  as.surface(gridList_big, surf_ns_k2add),
  col = plasma(256), 
  main = "",
  xlab = "Longitude",
  ylab = "Latitude",
  cex.lab = 1.2, 
  horizontal = TRUE,
  legend.shrink = 0.7, 
  legend.width = 0.8, 
  legend.args = list(
    text = expression(mu*g/m^{3}),
    side = 4,  las = 1, line = 1, cex = 1.0
  )
  
)
world(add = TRUE, col = "black", lwd = 1)
points(
  x = stations_df_day$Longitude,
  y = stations_df_day$Latitude,
  pch = 20,        # "x"
  col = "black",
  cex = 0.2,
  lwd = 0.5
)

library(cmocean)
imagePlot(
  as.surface(gridList_big, csim_ns_k2add$SE),
  col = mako(256),
  main = "",
  xlab = "Longitude",
  ylab = "", 
  yaxt = "n",
  cex.lab = 1.2,
  horizontal = TRUE,
  legend.shrink = 0.7,
  legend.width = 0.8, 
  legend.args = list(
    text = expression(mu*g/m^{3}),
    side = 4,  las = 1, line = 1, cex = 1.0
  )
)
world(add = TRUE, col = "black", lwd = 1)
points(
  x = stations_df_day$Longitude,
  y = stations_df_day$Latitude,
  pch = 16,        # "x"
  col = "black",
  cex = 0.2,
  lwd = 0.5
)
dev.off()
