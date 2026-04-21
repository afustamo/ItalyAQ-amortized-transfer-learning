# importing parameters from STUN ####
library(rhdf5)
library(LatticeKrig)
library(tictoc)
setwd("C:/Users/anton/Desktop/Research/Italy_AQ_paper")

#import hacked LK files for additional kappa2 estimation
source("Italy_AQ_AmortisedLatticeKrig/R_scripts/hacked_LatticeKrig.R")

AQ_file_path <- "Italy_AQ_AmortisedLatticeKrig/data/IAQ_NO2_resid_2019.h5"
df_AQ <- h5read(AQ_file_path, "fields")
lat <- rev(h5read(AQ_file_path, "lat"))
lon <- rev(h5read(AQ_file_path, "lon"))
time <- h5read(AQ_file_path, "time")

# load in i2i network clim outputs
file_path <- "Italy_AQ_AmortisedLatticeKrig/results/italy_aq_outputs/df_rscale_30rep.h5"
h5ls(file_path)

# parameters from global stun estimator
df_STUN <- h5read(file_path, "stun")

# extract/make params for STUN
kappa2_s <- exp(df_STUN[,,1,])
awght_s <- kappa2_s + 4
# need to transform for LK 
theta_s <- -df_STUN[,,2,] + pi/2
# for plotting intuitive angles
# theta_s <- -df_STUN[,,2]
rho_s <- df_STUN[,,3,]

# shift everything to look nice on plots
df_AQ <- aperm(df_AQ, c(2,1,3))
df_AQ <- df_AQ[, dim(df_AQ)[2]:1, ]
kappa2_s <- kappa2_s[, dim(kappa2_s)[2]:1, ]
theta_s <- theta_s[, dim(theta_s)[2]:1, ]
rho_s <- rho_s[, dim(rho_s)[2]:1, ]

# for later coding into LK
rhox_s <- sqrt(rho_s)
rhoy_s <- 1/rhox_s

# okay now here's how we use these parameters in LK:
# let's just use the STUN ones for now 
# and let's use the first day 
time_sel <- 1
kappa2_1 <- kappa2_s[,,time_sel]
theta_1 <- theta_s[,,time_sel]
rho_1 <- rho_s[,,time_sel]
rhox_1 <- sqrt(rho_1)
rhoy_1 <- 1/rhox_1

# some grid setup 6.2 35.1
# [2,] 18.9 47.8
# rows <- dim(df_AQ)[1]
rows_ctm <- length(seq( 6.2,18.9,by=.1))
gridList_ctm<- list( x= seq( 6.2,18.9,by=.1),
                 y= seq( 35.1,47.8,by=.1) )
sGrid_ctm<- make.surface.grid(gridList_ctm)

# create H tensor out of params
H11 <- ( rhox_1^2 * (cos(theta_1))^2) + ( rhoy_1^2 * (sin(theta_1))^2 ) 
H12 <- (rhoy_1^2 - rhox_1^2)*(sin(theta_1)*cos(theta_1))
H21 <- H12 
H22 <- (rhox_1^2 * (sin(theta_1))^2) + (rhoy_1^2 * (cos(theta_1))^2)

# fill the high dimensional stencil (9 fields)
stencil_tensor <- array( NA, c( rows_ctm,rows_ctm,9))
stencil_tensor[,,1] <- 0.5 * H12
stencil_tensor[,,2] <- -H22
stencil_tensor[,,3] <- -0.5 * H12
stencil_tensor[,,4] <- -H11
stencil_tensor[,,5] <- kappa2_1 + 2 * H11 + 2 * H22
stencil_tensor[,,6] <- -H11
stencil_tensor[,,7] <- -0.5 * H12
stencil_tensor[,,8] <- -H22
stencil_tensor[,,9] <- 0.5 * H12

# next, we put everything into awght obj of a particular class
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

awght_obj <- list( x= gridList_ctm$x,  y= gridList_ctm$y, z=stencil_tensor )

class( awght_obj)<- "multivariateSurfaceGrid"

# now inside of your LKinfo, set this to be your awght object, like so:
LKinfo_ns <- LKrigSetup(sGrid_ctm, NC =rows_ctm, #rows
                     nlevel = 1,
                     a.wghtObject =  awght_obj,
                     normalize=FALSE,
                     NC.buffer = 0, overlap = 2.5, nu = 1)

LKinfo_norm <- LKrigSetup(sGrid_ctm, NC =rows_ctm, #rows
                     nlevel = 1,
                     a.wghtObject =  awght_obj,
                     normalize=TRUE,
                     NC.buffer = 0, overlap = 2.5, nu = 1)

nonstat_Q <- LKrig.precision(LKinfo_norm)

# estimating Lattice ####
## importing data ####

# d <- as.Date(time)[time_sel]+1
d <- as.Date("2019-01-01")
HPC <- F
setwd("C:/Users/anton/Desktop/Research/Italy_AQ_paper/Archivio")

d <- d + 1
p <- "AQ-FRK/v.3.0.6/"
source(paste0(p, "Z_settings/Z_settings.R"), local = TRUE)

p <- "AQ-FRK/v.3.0.1/"
source(paste0(p, "A_input/A_input.R"), local = TRUE)

p <- "AQ-FRK/v.3.0.1/"
source(paste0(p, "C_BAUs/C_BAUs.R"), local = TRUE)

# adding lag data
d <- d - 1
p <- "AQ-FRK/v.3.0.6/"
source(paste0(p, "Z_settings/Z_settings.R"), local = TRUE)

p <- "AQ-FRK/v.3.0.1/"
crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
EEA_t1 <- F
if (EEA_t1) {
  if (!HPC) {
    p <- "AQ-FRK/v.3.0.0/"
    load(paste0(p, "A_input/data/input/AQ_EEA_v100_ST.rda"))
    slot(AQ_EEA_v100_ST@sp, "proj4string") <- crs_wgs84
    A1_aq_eea_t1 <- #temporal subset
      AQ_EEA_v100_ST[, which(index(AQ_EEA_v100_ST@time) %in% sub_t[2])] 
    A1_aq_eea_t1 <- A1_aq_eea_t1[A1_aq_eea_t1@coords[, 1] > bbox[1] + .25 &
                                   A1_aq_eea_t1@coords[, 1] < bbox[2] - .25 &
                                   A1_aq_eea_t1@coords[, 2] > bbox[3] + .25 &
                                   A1_aq_eea_t1@coords[, 2] < bbox[4] - .25, ]
  }
}

if (!HPC) {
  p <- "AQ-FRK/v.3.0.0/"
  load(paste0(p, "A_input/data/input/AQ_CAMS_v100_ST.rda"))
}
A2_aq_cams_t1 <- #temporal subset
  AQ_CAMS_v100_ST[, which(index(AQ_CAMS_v100_ST@time) %in% sub_t[2])] #v.3.0.1
#spatial subset
A2_aq_cams_t1 <- A2_aq_cams_t1[A2_aq_cams_t1@coords[, 1] > bbox[1] &
                                 #maybe change?
                                 A2_aq_cams_t1@coords[, 1] < bbox[2] &
                                 A2_aq_cams_t1@coords[, 2] > bbox[3] &
                                 A2_aq_cams_t1@coords[, 2] < bbox[4], ]
names(A2_aq_cams_t1@data)[4] <-
  "AQ_CAMS_NO2"
rm(AQ_CAMS_v100_ST,AQ_EEA_v100_ST)
#
p <- "AQ-FRK/v.3.0.1/"

setwd("C:/Users/anton/Desktop/Research/Italy_AQ_paper/Italy_AQ_AmortisedLatticeKrig")

# save.image(file = paste0("data/images/input_data_preLattice_",d,".RData"))

## preparing input datasets ####
# load(file = paste0("data/images/input_data_preLattice_",d,".RData"))
Z <- as.matrix(over(A1_aq_eea, C2_BAUs_ext)[, c(5:7, 10, 11, 12, 14)])
Z_a2 <- as.matrix(over(A1_aq_eea, A2_aq_cams[, 4]))
Z <- cbind(Z, Z_a2)
colnames(Z)[8] <- "AQ_CAMS_NO2_t1"
# Z <- scale(Z)
# or leave unsclaed , dividing ssr by 10^6
Z[,"A3_ssr"] <- Z[,"A3_ssr"] / 10^6
head(Z)

if(EEA_t1){
Z <- cbind(Z, A1_aq_eea_t1@data$AQ_EEA_NO2)
colnames(Z)[9] <- "eea_t1"}

A1_aq_eea_df <- A1_aq_eea@data
s <- A1_aq_eea_df[, c(1, 2)]
y <-  A1_aq_eea_df$AQ_EEA_NO2

set.seed(as.numeric(sub_t[1]))

test <- sample(1:length(y), length(y) / 10)
train <- 1:length(y)
train <- train[-which(train %in% test)]
y_train <- y[train]
test <- test[order(test)]
y_test <- y[test]
s_train <- s[train, ]
s_test <- s[test, ]
Z_train <- Z[train, ]
Z_test <- Z[test, ]

if(EEA_t1){
naidx_test <- is.na(Z_test[, "eea_t1"])
y_test <- y_test[!naidx_test]
s_test <- s_test[!naidx_test, ]
Z_test <- Z_test[!naidx_test, ]
}

y_train_mod <- y_train
naidx <- is.na(y_train_mod)
y_train_mod <- y_train_mod[!naidx]
s_train_mod <- s_train[!naidx, ]
Z_train_mod <- Z_train[!naidx, ]

if(EEA_t1){
naidx <- is.na(Z_train_mod[, "eea_t1"])
y_train_mod <- y_train_mod[!naidx]
s_train_mod <- s_train_mod[!naidx, ]
Z_train_mod <- Z_train_mod[!naidx, ]
}


## traning models ####
tic()
Lmodel_0 <- LatticeKrig(s_train_mod, y_train_mod, Z = Z_train_mod)
toc()
Lmodel_0

tic()
Lmodel_0_norm <- LatticeKrig(s_train_mod, y_train_mod, Z = Z_train_mod, normalize = TRUE)
toc()
Lmodel_0_norm

# load(paste0("data/LKinfo_ctms_",d,".rda")) # d-1 for the image
tic()
Lmodel_ctms <- LatticeKrig(s_train_mod, y_train_mod, Z = Z_train_mod, 
                           LKinfo = LKinfo_ns) #used LKinfo_test the last time
toc()
Lmodel_ctms

tic()
Lmodel_ctms_norm <- LatticeKrig(s_train_mod, y_train_mod, Z = Z_train_mod, 
                                LKinfo = LKinfo_norm) #used LKinfo_test the last time
toc()
Lmodel_ctms_norm

tic()
LKinfo_ns_new <- LKinfo_ns
LKinfo_ns_new$a.wght <- NA
LKinfo_ns_new$a.wghtObject <- NULL
LKinfo_ns_new$a.wght <- list()
LKinfo_ns_new$a.wght[[1]] <- 4.01
attr(LKinfo_ns_new$a.wght, "isotropic") <- TRUE


# Lmodel_ctms_norm_add <- LatticeKrig_hacked(
#   s_train_mod,
#   y_train_mod,
#   Z = Z_train_mod,
#   nlevel = 1, 
#   NC = rows_ctm, 
#   normalize = TRUE, 
#   NC.buffer = 0, 
#   overlap = 2.5, 
#   findAwght = T
# )

Lmodel_ctms_norm_add <- LatticeKrig_hacked(
  s_train_mod,
  y_train_mod,
  Z = Z_train_mod,
  LKinfo = LKinfo_ns_new, 
  findAwght = T
)

LKinfo_norm_hack <- LKinfo_norm
LKinfo_norm_hack$a.wght[[1]][,5] <- LKinfo_norm_hack$a.wght[[1]][,5] + 
  (Lmodel_ctms_norm_add$LKinfo$a.wght[[1]]-4)

Lmodel_ctms_norm_hack <- LatticeKrig(s_train_mod, y_train_mod, Z = Z_train_mod, 
                                     LKinfo = LKinfo_norm_hack)
Lmodel_ctms_norm_hack
toc()





ncams <- which(colnames(Z_train_mod)=="AQ_CAMS_NO2_t1")
Lmodel_ctms_norm_add_nc <- LatticeKrig_hacked(
  s_train_mod,
  y_train_mod,
  Z = Z_train_mod[,-ncams],
  LKinfo = LKinfo_ns_new, 
  findAwght = T
)
LKinfo_norm_hack_nc <- LKinfo_norm
LKinfo_norm_hack_nc$a.wght[[1]][,5] <- LKinfo_norm_hack_nc$a.wght[[1]][,5] + 
  (Lmodel_ctms_norm_add_nc$LKinfo$a.wght[[1]]-4)

Lmodel_ctms_norm_hack_nc <- LatticeKrig(s_train_mod, y_train_mod, 
                                        Z = Z_train_mod[,-ncams], 
                                     LKinfo = LKinfo_norm_hack_nc)
Lmodel_ctms_norm_hack_nc





tic()
Lmodel_1 <- LatticeKrig(
  s_train_mod,
  y_train_mod,
  Z = Z_train_mod,
  nlevel = 1,
  NC = 128,
  normalize = FALSE,
  NC.buffer = 0,
  overlap = 2.5,
  findAwght = T
)
toc()
Lmodel_1

tic()
Lmodel_1_norm <- LatticeKrig(
  s_train_mod,
  y_train_mod,
  Z = Z_train_mod,
  nlevel = 1,
  NC = 128,
  normalize = TRUE,
  NC.buffer = 0,
  overlap = 2.5,
  findAwght = T
)
toc()
Lmodel_1_norm


## results ####
# coefficients
coef_t <- cbind(
  Lmodel_0$d.coef, Lmodel_0_norm$d.coef, 
  Lmodel_ctms$d.coef,Lmodel_ctms_norm$d.coef,
  Lmodel_1$d.coef, Lmodel_1_norm$d.coef
)
colnames(coef_t) <- c(
  "Lmodel_0","Lmodel_0_norm", 
  "Lmodel_ctms","Lmodel_ctms_norm",
  "Lmodel_1","Lmodel_1_norm"
)
coef_t
# 

# Log ligkelihood
logL_t <- cbind(
  Lmodel_0$lnProfileLike, Lmodel_0_norm$lnProfileLike,
  Lmodel_ctms$lnProfileLike, Lmodel_ctms_norm$lnProfileLike,
  Lmodel_1$lnProfileLike, Lmodel_1_norm$lnProfileLike
)
colnames(logL_t) <- c(
  "Lmodel_0","Lmodel_0_norm", 
  "Lmodel_ctms","Lmodel_ctms_norm",
  "Lmodel_1","Lmodel_1_norm"
)
logL_t

# effective degree of freedom
eff_df_t <- cbind(
  Lmodel_0$eff.df, Lmodel_0_norm$eff.df,
  Lmodel_ctms$eff.df, Lmodel_ctms_norm$eff.df,
  Lmodel_1$eff.df, Lmodel_1_norm$eff.df
)
colnames(eff_df_t) <- c(
  "Lmodel_0","Lmodel_0_norm", 
  "Lmodel_ctms","Lmodel_ctms_norm",
  "Lmodel_1","Lmodel_1_norm"
)
eff_df_t

# Tau
tau_t <- cbind(
  Lmodel_0$tau.MLE, Lmodel_0_norm$tau.MLE,
  Lmodel_ctms$tau.MLE, Lmodel_ctms_norm$tau.MLE,
  Lmodel_1$tau.MLE, Lmodel_1_norm$tau.MLE
)
colnames(tau_t) <- c(
  "Lmodel_0","Lmodel_0_norm", 
  "Lmodel_ctms","Lmodel_ctms_norm",
  "Lmodel_1","Lmodel_1_norm"
)
print(paste("std.dev of the process:",sd(y_train_mod)))
tau_t

#sigma
sigma_t <- cbind(
  Lmodel_0$sigma2.MLE, Lmodel_0_norm$sigma2.MLE,
  Lmodel_ctms$sigma2.MLE, Lmodel_ctms_norm$sigma2.MLE,
  Lmodel_1$sigma2.MLE, Lmodel_1_norm$sigma2.MLE
)
colnames(sigma_t) <- c(
  "Lmodel_0","Lmodel_0_norm", 
  "Lmodel_ctms","Lmodel_ctms_norm",
  "Lmodel_1","Lmodel_1_norm"
)
print(paste("std.dev of the process:",sd(y_train_mod)))
sigma_t

# summary(Lmodel_1$residuals)
# summary(Lmodel_ctms$residuals)
summary(y_train_mod)

summary(abs(Lmodel_0$residuals))
summary(abs(Lmodel_0_norm$residuals))
summary(abs(Lmodel_1$residuals))
summary(abs(Lmodel_1_norm$residuals))
summary(abs(Lmodel_ctms$residuals))
summary(abs(Lmodel_ctms_norm$residuals))
summary(abs(Lmodel_ctms_norm_hack$residuals))


# Surfaces ####

# reducing the area considered, because standard error is very slow!
C2_BAUs_ext_cams1 <- cbind(C2_BAUs_ext@data,over(C2_BAUs_ext,A2_aq_cams_t1))
Z_map <- C2_BAUs_ext_cams1[, c(5:7, 10, 11, 12, 14,19)]
row.names(Z_map)<-NULL


# Z_map <- scale(Z_map)

colnames(Z_map)[8] <- "AQ_CAMS_NO2_t1"
# Z_center <- attr(Z, "scaled:center")
# Z_scale  <- attr(Z, "scaled:scale")
# Z_map <- scale(Z_map, center = Z_center, scale = Z_scale) #ask?
Z_map[,"A3_ssr"] <- Z_map[,"A3_ssr"]/10^6
summary(Z_map)
grid <- C2_BAUs_ext@coords
row.names(grid)<-NULL

Z_map <- as.matrix(Z_map)

# surf_0 <- predict(Lmodel_0,grid,Z_map)

# surf_ctms <- predict(Lmodel_ctms,grid,Z_map)
x_c <- unique(grid[,1])
y_c <- unique(grid[,2])
gridList <- list(x=x_c,y=y_c)
gridxy <- make.surface.grid(gridList)


tic()
surf_1_full <- predict(Lmodel_1,gridxy,Z_map)
toc()
surf_1_mu <- predict(Lmodel_1,gridxy,Znew=Z_map,just.fixed=T)
surf_1_GP <- surf_1_full - surf_1_mu

surf_1_norm_full <- predict(Lmodel_1_norm,gridxy,Z_map)
surf_1_norm_mu <- predict(Lmodel_1_norm,gridxy,Znew=Z_map,just.fixed=T)
surf_1_norm_GP <- surf_1_norm_full - surf_1_norm_mu

surf_0_full <- predict(Lmodel_0,gridxy,Z_map)
surf_0_mu <- predict(Lmodel_0,gridxy,Znew=Z_map,just.fixed=T)
surf_0_GP <- surf_0_full - surf_0_mu

surf_0_norm_full <- predict(Lmodel_0_norm,gridxy,Z_map)
surf_0_norm_mu <- predict(Lmodel_0_norm,gridxy,Znew=Z_map,just.fixed=T)
surf_0_norm_GP <- surf_0_norm_full - surf_0_norm_mu


surf_ctms_full <- predict(Lmodel_ctms,gridxy,Z_map)
surf_ctms_mu <- predict(Lmodel_ctms,gridxy,Znew=Z_map,just.fixed=T)
surf_ctms_GP <- surf_ctms_full - surf_ctms_mu

surf_ctms_norm_full <- predict(Lmodel_ctms_norm,gridxy,Z_map)
surf_ctms_norm_mu <- predict(Lmodel_ctms_norm,gridxy,Znew=Z_map,just.fixed=T)
surf_ctms_norm_GP <- surf_ctms_norm_full - surf_ctms_norm_mu

surf_ctms_norm_hack_full <- predict(Lmodel_ctms_norm_hack,gridxy,Z_map)
surf_ctms_norm_hack_mu <- predict(Lmodel_ctms_norm_hack,gridxy,Znew=Z_map,just.fixed=T)
surf_ctms_norm_hack_GP <- surf_ctms_norm_hack_full - surf_ctms_norm_hack_mu 

surf_ctms_norm_hack_nc_full <- predict(Lmodel_ctms_norm_hack_nc,gridxy,Z_map[,-ncams])
surf_ctms_norm_hack_nc_mu <- predict(Lmodel_ctms_norm_hack_nc,gridxy,Znew=Z_map[,-ncams],just.fixed=T)
surf_ctms_norm_hack_nc_GP <- surf_ctms_norm_hack_nc_full - surf_ctms_norm_hack_nc_mu 


par(mfrow = c(2,3))
# make margins smaller
par(mar=c(2,2,2,2))
imagePlot(as.surface(gridList, surf_ctms_norm_full), axes = F, ylab="",
          xlab="Nonstationary (normalized)",
          zlim=c(0,50),col=plasma(16))
world(add = T)
imagePlot(as.surface(gridList, surf_ctms_norm_GP), axes = F, xlab="",ylab="",
          col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_ctms_norm_mu), axes = F, xlab="",ylab="",
          zlim=c(0,85),col=turbo(256))
world(add = T)



imagePlot(as.surface(gridList, surf_ctms_norm_hack_full), axes = F, ylab="",
          xlab="Nonstationary",
          zlim=c(0,85),col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_ctms_norm_hack_GP), axes = F, xlab="",ylab="",
          col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_ctms_norm_hack_mu), axes = F, xlab="",ylab="",
          zlim=c(0,85),col=turbo(256)) 
world(add = T)



imagePlot(as.surface(gridList, surf_ctms_norm_hack_nc_full), axes = F, ylab="",
          xlab="Nonstationary",
          zlim=c(0,85),col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_ctms_norm_hack_nc_GP), axes = F, xlab="",ylab="",
          col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_ctms_norm_hack_nc_mu), axes = F, xlab="",ylab="",
          zlim=c(0,85),col=turbo(256)) 
world(add = T)

par(mfrow = c(1,1))





zmin <- min(c(surf_ctms_norm_hack_full, surf_ctms_norm_hack_nc_full))
zmax <- max(c(surf_ctms_norm_hack_full, surf_ctms_norm_hack_nc_full))



par(mfrow = c(1,2))
imagePlot(as.surface(gridList, surf_ctms_norm_hack_full), axes = F, ylab="",
          xlab="with cams",
          zlim=c(zmin,zmax),col=plasma(16))
# world(add = T)
imagePlot(as.surface(gridList, surf_ctms_norm_hack_nc_full), axes = F, ylab="",
          xlab="without cams",
          zlim=c(zmin,zmax),col=plasma(16))
# world(add = T)
par(mfrow = c(1,1))


# default 
par(mfrow = c(1,3))
imagePlot(as.surface(gridList, surf_0_full), axes = F, ylab="",
          xlab="LK default mean", zlim=c(0,85),col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_0_GP), axes = F, ylab="",
          xlab="LK default GP",col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_0_mu), axes = F, ylab="",
          xlab="LK default prediction ",
          zlim=c(0,85),col=turbo(256))
world(add = T)
par(mfrow = c(1,1))



# stat and nonstat comparison
par(mfrow = c(4,3))
imagePlot(as.surface(gridList, surf_1_full), axes = F, ylab="",
          xlab="Stationary Mean",
          zlim=c(0,85),col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_1_GP), axes = F, ylab="",
          xlab="Stationary GP",
          col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_1_mu), axes = F, ylab="",
          xlab="Stationary Prediction",
          zlim=c(0,85),col=turbo(256)) 
world(add = T)


imagePlot(as.surface(gridList, surf_1_norm_full), axes = F, ylab="",
          xlab="Stationary (normalized)",
          zlim=c(0,85),col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_1_norm_GP), axes = F, xlab="",ylab="",
          col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_1_norm_mu), axes = F, xlab="",ylab="",
          zlim=c(0,85),col=turbo(256)) 
world(add = T)



imagePlot(as.surface(gridList, surf_ctms_full), axes = F, ylab="",
          xlab="Nonstationary",
          zlim=c(0,85),col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_ctms_GP), axes = F, xlab="",ylab="",
          col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_ctms_mu), axes = F, xlab="",ylab="",
          zlim=c(0,85),col=turbo(256)) 
world(add = T)

imagePlot(as.surface(gridList, surf_ctms_norm_full), axes = F, ylab="",
          xlab="Nonstationary (normalized)",
          zlim=c(0,85),col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_ctms_norm_GP), axes = F, xlab="",ylab="",
          col=turbo(256))
world(add = T)
imagePlot(as.surface(gridList, surf_ctms_norm_mu), axes = F, xlab="",ylab="",
          zlim=c(0,85),col=turbo(256)) 
world(add = T)

par(mfrow = c(1,1))



# okay now let's look at the differences 
par(mfrow = c(2,3))
imagePlot(as.surface(gridList,surf_1_full),axes=F,xlab="Stat",ylab="",
          col=turbo(256))
world(add=T)
imagePlot(as.surface(gridList,surf_ctms_full),axes=F,xlab="Nonstat",ylab="",
          col=turbo(256))
world(add=T)
image.plot(as.surface(gridList,surf_ctms_full-surf_1_full),
           axes=F,xlab="Difference",ylab="",col=turbo(256))
world(add=T)

imagePlot(as.surface(gridList,surf_1_norm_full),axes=F,xlab="Stat (Norm)",ylab="",
          col=turbo(256))
world(add=T)
imagePlot(as.surface(gridList,surf_ctms_norm_full),axes=F,xlab="Nonstat (Norm)",ylab="",
          col=turbo(256))
world(add=T)
image.plot(as.surface(gridList,surf_ctms_norm_full-surf_1_norm_full),
           axes=F,xlab="Difference (Norm)",ylab="",col=turbo(256))
world(add=T)
par(mfrow = c(1,1))


#what are the spatial process absolute values? 
par(mfrow = c(2,2))
imagePlot(as.surface(gridList,abs(surf_1_GP)),axes=F,xlab="Abs(Stat GP)",ylab="",
          col=turbo(256))
world(add=T, col = "white", lwd = 1.5)
imagePlot(as.surface(gridList,abs(surf_ctms_GP)),axes=F,xlab="Abs(Nonstat GP)",ylab="",
          col=turbo(256))
world(add=T, col = "white", lwd = 1.5)

imagePlot(as.surface(gridList, abs(surf_1_norm_GP)),axes=F, xlab="Abs(Stat GP) (Norm)",ylab="",
          col=turbo(256))
world(add = T, col = "white", lwd = 1.5)

imagePlot(as.surface(gridList, abs(surf_ctms_norm_GP)),axes=F, xlab="Abs(Nonstat GP) (Norm)",ylab="",
          col=turbo(256))
world(add = T, col = "white", lwd = 1.5)
par(mfrow = c(1,1))



tic()
csim_1 <- LKrig.sim.conditional(Lmodel_1, M = 100, x.grid = gridxy, Z.grid = Z_map)
toc()
tic()
csim_1_norm <- LKrig.sim.conditional(Lmodel_1_norm, M = 100, x.grid = gridxy, Z.grid = Z_map)
toc()
tic()
csim_ctms <- LKrig.sim.conditional(Lmodel_ctms, M = 100, x.grid = gridxy, Z.grid = Z_map)
toc()
tic()
csim_ctms_norm <- LKrig.sim.conditional(Lmodel_ctms_norm, M = 100, x.grid = gridxy, Z.grid = Z_map)
toc()



plotcol <- turbo(256)
par(mfrow = c(4,3))
image.plot(as.surface(gridxy, csim_1$ghat), col = plotcol, xlab = "Stat Csim Mean")
image.plot(as.surface(gridxy, csim_1$g.draw[,77]), col = plotcol, xlab = "Stat Sample Draw")
image.plot(as.surface(gridxy, csim_1$SE), col = plasma(256), xlab = "Stat SE")

image.plot(as.surface(gridxy, csim_1_norm$ghat), col = plotcol, xlab = "Stat Norm")
image.plot(as.surface(gridxy, csim_1_norm$g.draw[,77]), col = plotcol)
image.plot(as.surface(gridxy, csim_1_norm$SE), col = plasma(256))

image.plot(as.surface(gridxy, csim_ctms$ghat), col = plotcol, xlab = "Nonstat")
image.plot(as.surface(gridxy, csim_ctms$g.draw[,77]), col = plotcol)
image.plot(as.surface(gridxy, csim_ctms$SE), col = plasma(256))

image.plot(as.surface(gridxy, csim_ctms_norm$ghat), col = plotcol, xlab = "Nonstat Norm")
image.plot(as.surface(gridxy, csim_ctms_norm$g.draw[,77]), col = plotcol)
image.plot(as.surface(gridxy, csim_ctms_norm$SE), col = plasma(256))
par(mfrow = c(1,1))
