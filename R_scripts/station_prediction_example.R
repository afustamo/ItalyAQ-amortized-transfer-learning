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
library(tictoc)
library(scico)
library(here)
library(cmocean)
library(maps)

# ------------------------------------
# Choose which day
# ------------------------------------
chosen_day <- 305 #14 #305 is nov 1. 103 is apr 13

# import helpful funcs
source(here("R_scripts", "helpful_functions.R"))

set.seed(66)

# create grid 
rows_small <- length(seq(6.1, 18.8, by = .1))

gridList_small <- list(
  x = seq(6.1, 18.8, by = .1),
  y = seq(35.2, 47.9, by = .1)
)
sGrid_small <- make.surface.grid(gridList_small)


# ------------------------------------
# Dataset setup 
# ------------------------------------
stations_df <- readRDS(here("data", "eea_df_2023.rds"))
stations_df[, "ssr"] <- stations_df[, "ssr"] / 10^6

# for weights 
cams_df <- readRDS(here("data", "cams_df_2023padded.rds"))
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
  file_path   = here("results", "STUN_param_df_2023padded.h5"),
  dataset_name= "arx1_surround_30rep_output",
  day         = chosen_day,
  gridlist    = gridList_small,
  normalize   = TRUE,
  sanity_plotting = FALSE,
  sanity_sim      = FALSE
)
nonstat_Q <- LKrig.precision(LKinfo_ns)


# nonstat model with additonal kappa2 estimated 
# default
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
preds_df <- readRDS(here("data", "pred_df_2023.rds"))
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

summary(csim_stat$SE)
summary(csim_ns_k2add$SE)
summary(csim_ns$SE)


pdf(here("figures", "fine_grid_pred_nov1.pdf"), width = 7, height = 5)
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
