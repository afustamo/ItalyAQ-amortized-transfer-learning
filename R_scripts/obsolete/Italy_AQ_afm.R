stplot_y <- function(y,name,name_plot="NO2",namefile="plot_30days/plot.pdf"){
  Y <- y
  df <- data.frame(
    Longitude = rep(rep(dimnames(Y)[[2]], each = dim(Y)[1]), dim(Y)[3]),
    Latitude = rep(dimnames(Y)[[1]], dim(Y)[2] * dim(Y)[3]),
    time = rep(as.Date(as.numeric(dimnames(
      Y
    )[[3]])), each = dim(Y)[2] * dim(Y)[1]),
    var_i = c(Y)
  )
  names(df)[4]<-name
  df$Longitude <- as.numeric(df$Longitude)
  df$Latitude <- as.numeric(df$Latitude)
  centre <- unique(df[,c(1,2)])
  coordinates(centre) <- c("Longitude", "Latitude")
  gridded(centre) <- TRUE
  colnames(centre@coords) <- c("coords.x1", "coords.x2")
  crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
  slot(centre, "proj4string") <- crs_wgs84
  st <- STFDF(sp = centre,
                           time = unique(df$time),
                           data = df)
  pdf(file = namefile,width = 8,height = 8)
  print(stplot(st[,,name],main=name_plot))
  dev.off()
}

spplot_y <- function(y,name,name_plot="NO2",namefile="plot_30days/plot.pdf"){
  Y <- y
  df <- data.frame(
    Longitude = rep(dimnames(Y)[[2]], each = dim(Y)[1]),
    Latitude = rep(dimnames(Y)[[1]], dim(Y)[2]),
    var_i = c(Y)
  )
  names(df)[3]<-name
  df$Longitude <- as.numeric(df$Longitude)
  df$Latitude <- as.numeric(df$Latitude)
  centre <- unique(df[,c(1,2)])
  coordinates(centre) <- c("Longitude", "Latitude")
  gridded(centre) <- TRUE
  colnames(centre@coords) <- c("coords.x1", "coords.x2")
  crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
  slot(centre, "proj4string") <- crs_wgs84
  sp <- SpatialPixelsDataFrame(centre,df)
  pdf(file = namefile,width = 8,height = 8)
  print(spplot(sp[,name],main=name_plot))
  dev.off()
}

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/Italy_AQ_AmortisedLatticeKrig")
library(spacetime)
library(sp)
library(abind)

load("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-CAMS/v.1.0.0/data/output/no2_2023.rda")

stplot_y(y[,,1:30],"NO2","raw data","plot_30days/raw_data.pdf")

#norm1: compute difference between raw data and monthly average, and divide by std####
y_avg <- y
for (i in 1:15) {
  y_avg[,,i] <- apply(y[,,1:30],c(1,2),mean)
}
for (i in 16:45) {
  y_avg[,,i] <- apply(y[,,(i-15):(i+14)],c(1,2),mean)
}

y_sd <- y
for (i in 1:15) {
  y_sd[,,i] <- apply(y[,,1:30],c(1,2),sd)
}
for (i in 16:45) {
  y_sd[,,i] <- apply(y[,,(i-15):(i+14)],c(1,2),sd)
}

y_norm1 <- y
y_norm1[,,1:59] <- y[,,1:59] - y_avg[,,1:59]
y_norm1[,,1:59] <- y_norm1[,,1:59] / y_sd[,,1:59]

stplot_y(y_norm1[,,1:30],"NO2","norm1: normalise data","plot_30days/norm1.pdf")
stplot_y(y_avg[,,1:30],"NO2","norm1: avg_month","plot_30days/norm1_avg.pdf")
stplot_y(y_sd[,,1:30],"NO2","norm1: sd_month","plot_30days/norm1_sd.pdf")

#norm2: difference day by day ####
y_diff <- apply(y[,,1:59],c(1,2),diff)
dim(y_diff)
y_diff <- aperm(y_diff,perm = c(2,3,1))
dimnames(y_diff) <- dimnames(y[,,1:58])
stplot_y(y_diff[,,31:45],"NO2_diff","norm2: diff NO2","plot_30days/norm2_diff.pdf")
y_diff_sd <- y_diff
for (i in 1:15) {
  y_diff_sd[,,i] <- apply(y_diff[,,1:30],c(1,2),sd)
}
for (i in 16:44) {
  y_diff_sd[,,i] <- apply(y_diff[,,(i-15):(i+14)],c(1,2),sd)
}
stplot_y(y_diff_sd[,,1:30],"NO2_sd","norm2: sd diff NO2","plot_30days/norm2_sd.pdf")
y_norm2 <- y_diff
for (i in 1:58) {
  y_norm2[,,i] <- y_diff[,,i]/y_diff_sd[,,i]
}
stplot_y(y_norm2[,,31:45],"NO2","norm2: (raw[t] - raw[t-1]) / sd(raw[t]-raw[t-1])","plot_30days/norm2.pdf")

# norm3: linear trend weather + emissions####
load("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-CAMS/v.1.0.0/data/output/ST/AQ_CAMS_v100_ST_no2.rda")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-CAMS/v.1.0.0/data/output/no2_2023.rda")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-FRK/v.3.0.0/A_input/data/input/WE_C3S_v100_ST_ERA5SL.rda")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-FRK/v.3.0.0/A_input/data/input/WE_C3S_v100_ST_ERA5Land.rda")
EM_CAMS_v200_ST <- readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-FRK/v.3.0.8/A_input/data/input/EM_CAMS_v200_ST.rds")

t <- which(index(AQ_CAMS_v100_ST@time)=="2019/01/01")
stplot(AQ_CAMS_v100_ST[,1:3,"NO2"])
crs_wgs84 <- CRS(SRS_string = "EPSG:4326") #reference system for WGS 84
slot(AQ_CAMS_v100_ST@sp, "proj4string") <- crs_wgs84

stplot(WE_C3S_v100_ST_ERA5Land[,1:3,"t2m"])
slot(WE_C3S_v100_ST_ERA5Land@sp, "proj4string") <- crs_wgs84

stplot(WE_C3S_v100_ST_ERA5SL[,1:3,"t2m"])
slot(WE_C3S_v100_ST_ERA5SL@sp, "proj4string") <- crs_wgs84

stplot(EM_CAMS_v200_ST[,1:3,"EM_CAMS_NO2"])
EM_CAMS_v200_ST@data$EM_CAMS_NO2 <- scale(EM_CAMS_v200_ST@data$EM_CAMS_NO2)
stplot(EM_CAMS_v200_ST[,1:3,"EM_CAMS_NO2"])
slot(EM_CAMS_v200_ST@sp, "proj4string") <- crs_wgs84

cams_era5land <- over(AQ_CAMS_v100_ST[,t:(t+99)],WE_C3S_v100_ST_ERA5Land[,t:(t+99)])
cams_era5sl <- over(AQ_CAMS_v100_ST[,t:(t+99)],WE_C3S_v100_ST_ERA5SL[,t:(t+99)])
cams_em <- over(AQ_CAMS_v100_ST[,t:(t+99)],EM_CAMS_v200_ST[,1:100])

AQ_CAMS_v100_ST <- AQ_CAMS_v100_ST[,t:(t+99)]

cams_era5land <- cbind(AQ_CAMS_v100_ST@data,cams_era5land)
cams_era5sl <- cbind(AQ_CAMS_v100_ST@data,cams_era5sl)
cams_em <- cbind(AQ_CAMS_v100_ST@data,cams_em)
cams_em$EM_CAMS_NO2[cams_em$EM_CAMS_NO2>quantile(cams_em$EM_CAMS_NO2,probs = .95)] <- quantile(cams_em$EM_CAMS_NO2,probs = .95)
summary(cams_em)
names(cams_era5sl)<-paste0("sl_",names(cams_era5sl))
cams_we <- cbind(cams_era5land,cams_era5sl[,-c(1:7)])
for (i in 8:15) {
  var <- names(cams_we)[i]
  n_idx <- grep(paste0("sl_",var),names(cams_we))
  cams_we[is.na(cams_we[,i]),i] <- cams_we[is.na(cams_we[,i]),n_idx]
}

cams_we <- cams_we[,c(1:4,8:13,15:16)]
cams_we_em <- cbind(cams_we,cams_em[,-c(1:7)])
names(cams_we_em)[13]<-"EM_NO2"
cams_we_em$lag_cams_no2 <- c(rep(0,16900),cams_we_em$NO2[-c((1690000-16899):1690000)])

#skip to lm1
lm0 <- lm(NO2 ~ .-Longitude - Latitude - time,data=cams_we_em)
summary(lm0)
lm1 <- lm(NO2 ~ .-Longitude - Latitude - time - lag_cams_no2 ,data=cams_we_em)
summary(lm1)

cams_we_em$predict_lm <- lm0$fitted.values
cams_we_em$residuals_ls <- cams_we_em$NO2 - cams_we_em$predict_lm
cams_we_em$predict_lm1 <- lm1$fitted.values
cams_we_em$residuals_ls1 <- cams_we_em$NO2 - cams_we_em$predict_lm1

AQ_CAMS_v100_ST_mod <- AQ_CAMS_v100_ST
centre <- unique(cams_we[,c(1:2)])
centre[,1] <- as.numeric(centre[,1])
centre[,2] <- as.numeric(centre[,2])
coordinates(centre) <- c("Longitude", "Latitude")
gridded(centre) <- TRUE
colnames(centre@coords) <- c("coords.x1", "coords.x2")
crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
slot(centre, "proj4string") <- crs_wgs84
cams_we_st <- STFDF(sp = centre,
                         time = as.Date(unique(cams_we_em$time)),
                         data = cams_we_em)
for (i in sample(1:99,4)) {
  print(spplot(cams_we_st[,i,c("residuals_ls","residuals_ls1")],main="ARX1 residual 1-feb 15-feb 2019",at=seq(-30,50,length.out=100))
)
}
spplot(cams_we_st[,32,c("residuals_ls","residuals_ls1")],main="ARX1 residual 1-feb 15-feb 2019",at=seq(-30,50,length.out=100))

pdf("plot/residuals_ARX1.pdf")
stplot(cams_we_st[,32:47,"residuals_ls"],main="ARX1 residual 1-feb 15-feb 2019")
cams_we_st@data$residuals_ls_scale <- scale(cams_we_st@data$residuals_ls)
stplot(cams_we_st[,32:47,"residuals_ls_scale"],main="SCALED ARX1 residual 1-feb 15-feb 2019")
dev.off()

# 
AQ_CAMS_v100_ST_mod <- AQ_CAMS_v100_ST
centre <- unique(cams_we[,c(1:2)])
centre[,1] <- as.numeric(centre[,1])
centre[,2] <- as.numeric(centre[,2])
coordinates(centre) <- c("Longitude", "Latitude")
gridded(centre) <- TRUE
colnames(centre@coords) <- c("coords.x1", "coords.x2")
crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
slot(centre, "proj4string") <- crs_wgs84
cams_we_st <- STFDF(sp = centre,
                    time = as.Date(unique(cams_we_em$time)),
                    data = cams_we_em)
pdf("plot/residuals_ARX1.pdf")
stplot(cams_we_st[,32:47,"residuals_ls"],main="ARX1 residual 1-feb 15-feb 2019",at=seq(-30,30,length.out=100))
cams_we_st@data$residuals_ls_scale <- scale(cams_we_st@data$residuals_ls)
stplot(cams_we_st[,32:47,"residuals_ls_scale"],main="SCALED ARX1 residual 1-feb 15-feb 2019")
dev.off()

#
summary(cams_we_st@data$residuals_ls)
sd(lm0$residuals[1:16900])
hist(cams_we_st@data$residuals_ls)

saveRDS(cams_we_em,file = "data/cams_we_em_2019.rds")

# plot again
library(sp)
library(spacetime)
cams_we_em_2019 <- readRDS("data/cams_we_em_2019.rds")

centre <- unique(cams_we_em_2019[,c(1:2)])
centre[,1] <- as.numeric(centre[,1])
centre[,2] <- as.numeric(centre[,2])
coordinates(centre) <- c("Longitude", "Latitude")
gridded(centre) <- TRUE
colnames(centre@coords) <- c("coords.x1", "coords.x2")
crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
slot(centre, "proj4string") <- crs_wgs84
cams_we_st <- STFDF(sp = centre,
                    time = as.Date(unique(cams_we_em_2019$time)),
                    data = cams_we_em_2019)
stplot(cams_we_st[,1:3,"residuals_ls"])


