# importing parameters from STUN ####
library(rhdf5)

AQ_file_path <- "data/IAQ_NO2_resid_2019.h5"
df_AQ <- h5read(AQ_file_path, "fields")
lat <- rev(h5read(AQ_file_path, "lat"))
lon <- rev(h5read(AQ_file_path, "lon"))
time <- h5read(AQ_file_path, "time")

# load in i2i network clim outputs
file_path <- "data/df_rscale_30rep.h5"
h5ls(file_path)

# parameters from global stun estimator
df_STUN <- h5read(file_path, "stun")

# extract/make params for STUN
kappa2_s <- exp(df_STUN[,,1,])
awght_s <- kappa2_s + 4
# need to transform for LK 
theta_s <- pi/2 - df_STUN[,,2,]
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

# H11 <- 1
# H12 <- 0
# H21 <- 0
# H22 <- 1
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
# 
# stencil_tensor <- array( NA, c( rows,rows,9))
# stencil_tensor[,,1] <- 0
# stencil_tensor[,,2] <- -1
# stencil_tensor[,,3] <- 0
# stencil_tensor[,,4] <- -1
# # stencil_tensor[,,5] <- kappa2_1 + 2 * H11 + 2 * H22
# stencil_tensor[,,5] <- 4.317692
# stencil_tensor[,,6] <- -1
# stencil_tensor[,,7] <- 0
# stencil_tensor[,,8] <- -1
# stencil_tensor[,,9] <- 0

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
LKinfo <- LKrigSetup(sGrid_ctm, NC =rows_ctm, #rows
                     nlevel = 1, 
                     a.wghtObject =  awght_obj, 
                     normalize=FALSE, 
                     NC.buffer = 0, overlap = 2.5, nu = 1) 

save(LKinfo,file = paste0("data/LKinfo_ctms_",as.Date(time)[time_sel],".rda"))
save.image(file = paste0("data/images/LKinfo_",as.Date(time)[time_sel],".RData"))
rm(list = setdiff(ls(),c("time","time_sel")))
gc()

# estimating Lattice ####
## importing data ####

# d <- as.Date(time)[time_sel]+1
d <- as.Date("2019-01-01")
HPC <- F
setwd(
  "~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2"
)
d <- d+1
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
if(EEA_t1){
if (!HPC) {
  p <- "AQ-FRK/v.3.0.0/"
  load(paste0(p, "A_input/data/input/AQ_EEA_v100_ST.rda"))
}
slot(AQ_EEA_v100_ST@sp, "proj4string") <- crs_wgs84
A1_aq_eea_t1 <- #temporal subset
  AQ_EEA_v100_ST[, which(index(AQ_EEA_v100_ST@time) %in% sub_t[2])] 
A1_aq_eea_t1 <- A1_aq_eea_t1[A1_aq_eea_t1@coords[, 1] > bbox[1] + .25 &
                               A1_aq_eea_t1@coords[, 1] < bbox[2] - .25 &
                               A1_aq_eea_t1@coords[, 2] > bbox[3] + .25 &
                               A1_aq_eea_t1@coords[, 2] < bbox[4] - .25, ]
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

setwd(
  "~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/Italy_AQ_AmortisedLatticeKrig"
)

save.image(file = paste0("data/images/input_data_preLattice_",d,".RData"))

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

A1_aq_eea_df <-
  A1_aq_eea@data
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
Lmodel_0 <- LatticeKrig(s_train_mod, y_train_mod, Z = Z_train_mod)
Lmodel_0

load(paste0("data/LKinfo_ctms_",d,".rda")) # d-1 for the image
Lmodel_ctms <- LatticeKrig(s_train_mod, y_train_mod, Z = Z_train_mod, LKinfo = LKinfo) #used LKinfo_test the last time

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

save(Lmodel_0,Lmodel_1,Lmodel_ctms,file = paste0("data/model_trained_",d,".rda"))
rm(list = setdiff(ls(),c("d","Lmodel_0","Lmodel_1","Lmodel_ctms",
                         "C2_BAUs_ext","A2_aq_cams_t1","y_train_mod")))

## results ####
# coefficients
coef_t <- cbind(Lmodel_0$d.coef,Lmodel_ctms$d.coef,Lmodel_1$d.coef)
colnames(coef_t) <- c("Lmodel_0","Lmodel_ctms","Lmodel_1")
coef_t
# 

# Log ligkelihood
logL_t <- cbind(Lmodel_0$lnProfileLike,Lmodel_ctms$lnProfileLike,Lmodel_1$lnProfileLike)
colnames(logL_t) <- c("Lmodel_0","Lmodel_ctms","Lmodel_1")
logL_t

# effective degree of freedom
eff_df_t <- cbind(Lmodel_0$eff.df,Lmodel_ctms$eff.df,Lmodel_1$eff.df)
colnames(eff_df_t) <- c("Lmodel_0","Lmodel_ctms","Lmodel_1")
eff_df_t

# Tau
tau_t <- cbind(Lmodel_0$tau.MLE,Lmodel_ctms$tau.MLE,Lmodel_1$tau.MLE)
colnames(tau_t) <- c("Lmodel_0","Lmodel_ctms","Lmodel_1")
print(paste("std.dev of the process:",sd(y_train_mod)))
tau_t

#sigma
sigma_t <- cbind(Lmodel_0$sigma2.MLE,Lmodel_ctms$sigma2.MLE,Lmodel_1$sigma2.MLE)
colnames(sigma_t) <- c("Lmodel_0","Lmodel_ctms","Lmodel_1")
print(paste("std.dev of the process:",sd(y_train_mod)))
sigma_t

summary(Lmodel_1$residuals)
summary(Lmodel_ctms$residuals)
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

# surf_0 <- predict(Lmodel_0,grid,Z_map)

# surf_ctms <- predict(Lmodel_ctms,grid,Z_map)
x_c <- unique(grid[,1])
y_c <- unique(grid[,2])
gridList <- list(x=x_c,y=y_c)
gridxy <- make.surface.grid(gridList)

Z_map <- as.matrix(Z_map)

surf_1 <- predict(Lmodel_1,gridxy,Z_map)
surf_b <- predict(Lmodel_1,gridxy,drop.Z=T)
surf_a <- predict(Lmodel_1,gridxy,Znew=Z_map,just.fixed=T)

surf_ctms_1 <- predict(Lmodel_ctms,gridxy,Z_map)
surf_ctms_b <- predict(Lmodel_ctms,gridxy,drop.Z=T)
surf_ctms_a <- predict(Lmodel_ctms,gridxy,Znew=Z_map,just.fixed=T)

pdf(paste0("plot/comparing_LmodelCTM_Lmodel1_",d,".pdf"))
set.panel(2,3)
imagePlot(as.surface(gridList,surf_a),axes=F,xlab="",ylab="",
          zlim=c(0,85),col=turbo(256))
world(add=T)
surf_lonlat <- cbind(1,gridxy) %*% Lmodel_1$d.coef[1:3]
imagePlot(as.surface(gridList,surf_b-surf_lonlat))
          # ,axes=F,xlab="",ylab="",
          # zlim=c(-5,5),col=turbo(256))
world(add=T)
imagePlot(as.surface(gridList,surf_1),axes=F,xlab="",ylab="",
          zlim=c(0,85),col=turbo(256))
world(add=T)


# set.panel(1,3)
imagePlot(as.surface(gridList,surf_ctms_a),axes=F,xlab="",ylab="",
          zlim=c(0,85),col=turbo(256))
world(add=T)
surf_ctms_lonlat <- cbind(1,gridxy) %*% Lmodel_ctms$d.coef[1:3]
imagePlot(as.surface(gridList,surf_ctms_b-surf_ctms_lonlat))
# ,axes=F,xlab="",ylab="",
# zlim=c(-5,5),col=turbo(256))
world(add=T)
imagePlot(as.surface(gridList,surf_ctms_1),axes=F,xlab="",ylab="",
          zlim=c(0,85),col=turbo(256))
world(add=T)
dev.off()
# surfSE_0 <- predictSE(Lmodel_0,grid,Z_map)
# surfSE_ctms <- predictSE(Lmodel_ctms,grid,Z_map)
# surfSE_1 <- predictSE(Lmodel_1,grid,Z_map)

#very slow!

