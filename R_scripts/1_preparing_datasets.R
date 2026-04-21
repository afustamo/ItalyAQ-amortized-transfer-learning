library(sp)
library(spacetime)
library(FRK)
library(raster)
library(terra)
library(abind)

# import data and cut 2023 #### 
t <- seq.Date(as.Date("2022-12-10"),as.Date("2023-12-31"),"days")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-CAMS/v.1.0.0/data/output/ST/AQ_CAMS_v100_ST_no2.rda")
# AQ_CAMS_v100_ST_0 <- AQ_CAMS_v100_ST[,which(index(AQ_CAMS_v100_ST@time)== t[1]-1)]
AQ_CAMS_v100_ST <- AQ_CAMS_v100_ST[,which(index(AQ_CAMS_v100_ST@time)%in% t)]
#adding 2024
cams_df <- AQ_CAMS_v100_ST@data
cams_2024_no2 <- readRDS("data/ext_22_24/cams_2024_no2.rds")
head(cams_2024_no2)
cams_2024_no2$lon <- cams_2024_no2$lon + .05
cams_2024_no2$lat <- cams_2024_no2$lat + .05
cams_2024_no2 <- cams_2024_no2[,c("lon","lat","days","mean_no2")]
names(cams_2024_no2)<-c("Longitude","Latitude","time","NO2")
setdiff(unique(cams_df$Longitude),unique(cams_2024_no2$Longitude))
setdiff(unique(cams_df$Latitude),unique(cams_2024_no2$Latitude))
cams_df <- rbind(cams_df,cams_2024_no2)
cams_df$Longitude <- as.numeric(cams_df$Longitude)
cams_df$Latitude <- as.numeric(cams_df$Latitude)
cams_df$time <- as.Date(cams_df$time)
cams_df$NO2 <- as.numeric(cams_df$NO2)
cams_df <- cams_df[order(cams_df$time,cams_df$Latitude,cams_df$Longitude),]
cams_sp <- unique(cams_df[,c("Longitude","Latitude")])
coordinates(cams_sp)<-c("Longitude","Latitude")
cams_st <- STFDF(cams_sp,unique(cams_df$time),cams_df)
rm(AQ_CAMS_v100_ST,cams_2024_no2,cams_df,cams_sp)
gc()

#sl
load("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-FRK/v.3.0.0/A_input/data/input/WE_C3S_v100_ST_ERA5SL.rda")
WE_C3S_v100_ST_ERA5SL <- WE_C3S_v100_ST_ERA5SL[,which(index(WE_C3S_v100_ST_ERA5SL@time)%in% t)]
we_df_sl <- WE_C3S_v100_ST_ERA5SL@data
we_df_sl$Longitude <- as.numeric(we_df_sl$Longitude)
we_df_sl$Latitude <- as.numeric(we_df_sl$Latitude)

lf <- list.files("data/ext_22_24/we_2024/sl","rds")
ndata <- length(lf)
for (i in 1:ndata) {
  y <- readRDS(paste0("data/ext_22_24/we_2024/sl/", lf[i]))
  if (i == 1) {
    Y <- y
  } else{
    Y <- abind(Y, y, along = 4)
  }
}
varnames <- gsub(".rds","",lf)
dimnames(Y)[[4]] <- varnames
rm(y)
we_df_sl_24 <- data.frame(
  Longitude = rep(rep(dimnames(Y)[[2]], each = dim(Y)[1]), dim(Y)[3]),
  Latitude = rep(dimnames(Y)[[1]], dim(Y)[2] * dim(Y)[3]),
  time = rep(as.Date(as.numeric(dimnames(
    Y
  )[[3]])), each = dim(Y)[2] * dim(Y)[1])
)
we_df_sl_24 <- cbind(we_df_sl_24,
                       matrix(c(Y), ncol = dim(Y)[4]))
names(we_df_sl_24)[-c(1:3)] <- varnames
head(we_df_sl_24)
we_df_sl_24$Longitude <- as.numeric(we_df_sl_24$Longitude)
we_df_sl_24$Latitude <- as.numeric(we_df_sl_24$Latitude)
head(we_df_sl_24)
head(we_df_sl)
setdiff(we_df_sl$Longitude,we_df_sl_24$Longitude)
setdiff(we_df_sl$Latitude,we_df_sl_24$Latitude)
names(we_df_sl_24)[which(names(we_df_sl_24)=="winddirection")] <- "winddir"
we_df_sl <- rbind(we_df_sl,we_df_sl_24)
we_df_sl$time <- as.Date(we_df_sl$time)
we_df_sl <- we_df_sl[order(we_df_sl$time,we_df_sl$Latitude,we_df_sl$Longitude),]
we_sp_sl <- unique(we_df_sl[,c("Longitude","Latitude")])
coordinates(we_sp_sl)<-c("Longitude","Latitude")
gridded(we_sp_sl) <- TRUE
we_st_sl <- STFDF(we_sp_sl,unique(we_df_sl$time),we_df_sl)
stplot(we_st_sl[,1:2,"t2m"])
spplot(we_st_sl[,length(unique(we_df_sl$time)),"t2m"])
head(we_df_sl)

load("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-FRK/v.3.0.0/A_input/data/input/WE_C3S_v100_ST_ERA5Land.rda")
WE_C3S_v100_ST_ERA5Land <- WE_C3S_v100_ST_ERA5Land[,which(index(WE_C3S_v100_ST_ERA5Land@time)%in% t)]
we_df_land <- WE_C3S_v100_ST_ERA5Land@data
we_df_land$Longitude <- as.numeric(we_df_land$Longitude)
we_df_land$Latitude <- as.numeric(we_df_land$Latitude)

lf <- list.files("data/ext_22_24/we_2024/land","rds")
ndata <- length(lf)
for (i in 1:ndata) {
  y <- readRDS(paste0("data/ext_22_24/we_2024/land/", lf[i]))
  if (i == 1) {
    Y <- y
  } else{
    Y <- abind(Y, y, along = 4)
  }
}
varnames <- gsub(".rds","",lf)
dimnames(Y)[[4]] <- varnames
rm(y)
we_df_land_24 <- data.frame(
  Longitude = rep(rep(dimnames(Y)[[2]], each = dim(Y)[1]), dim(Y)[3]),
  Latitude = rep(dimnames(Y)[[1]], dim(Y)[2] * dim(Y)[3]),
  time = rep(as.Date(as.numeric(dimnames(
    Y
  )[[3]])), each = dim(Y)[2] * dim(Y)[1])
)
we_df_land_24 <- cbind(we_df_land_24,
                                 matrix(c(Y), ncol = dim(Y)[4]))
names(we_df_land_24)[-c(1:3)] <- varnames
head(we_df_land_24)
we_df_land_24$Longitude <- as.numeric(we_df_land_24$Longitude)
we_df_land_24$Latitude <- as.numeric(we_df_land_24$Latitude)
head(we_df_land_24)
head(we_df_land)
setdiff(we_df_land$Longitude,we_df_land_24$Longitude)
setdiff(we_df_land$Latitude,we_df_land_24$Latitude)
names(we_df_land_24)[10] <- "winddir"
we_df_land <- rbind(we_df_land,we_df_land_24)
we_df_land$time <- as.Date(we_df_land$time)
we_df_land <- we_df_land[order(we_df_land$time,we_df_land$Latitude,we_df_land$Longitude),]
we_sp_land <- unique(we_df_land[,c("Longitude","Latitude")])
coordinates(we_sp_land)<-c("Longitude","Latitude")
gridded(we_sp_land) <- TRUE
we_st_land <- STFDF(we_sp_land,unique(we_df_land$time),we_df_land)
stplot(we_st_land[,1:2,"t2m"])
spplot(we_st_land[,length(unique(we_df_sl$time)),"t2m"])
rm(WE_C3S_v100_ST_ERA5Land,WE_C3S_v100_ST_ERA5SL,we_df_land,we_df_land_24,
   we_df_sl,we_df_sl_24,we_sp_land,we_sp_sl)
gc()

# mnox_all
nox_all <- readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/EM-CAMS/v.3.0.0/data.nosync/EM_NOx-2.rds")
t2 <- c(t,seq.Date(t[length(t)]+1,t[length(t)]+31,by="days"))
nox_all <- nox_all[,,as.Date(as.numeric(dimnames(nox_all)[[3]]),origin=as.Date("1850-01-01")) %in% t2]
max_thr <- quantile(nox_all, .99) # da mettere nel preprocessing !
nox_all[nox_all > max_thr] <- max_thr

dimnames(nox_all)[[3]]<-as.Date(as.numeric(dimnames(nox_all)[[3]]),origin=as.Date("1850-01-01"))
Y <- nox_all
# EM_CAMS_v300_ST <- readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/EM-CAMS/v.3.0.0/data.nosync/ST/EM_NO2_CAMS_v300_ST_2023.rds")
# EM_CAMS_v300_ST_24 <- readRDS("data/ext_22_24/EM_CAMS_v300_ST_ALL.rds")
# EM_CAMS_v300_ST_24@data <- EM_CAMS_v300_ST_24@data[,c("Longitude","Latitude","time","EM_NOx")]
# EM_CAMS_v300_ST_24 <- EM_CAMS_v300_ST_24[,which(index(EM_CAMS_v300_ST_24@time)%in% t)]
gc()
nox_df <- data.frame(
  Longitude = rep(rep(dimnames(Y)[[2]], each = dim(Y)[1]), dim(Y)[3]),
  Latitude = rep(dimnames(Y)[[1]], dim(Y)[2] * dim(Y)[3]),
  time = rep(dimnames(
    Y
  )[[3]], each = dim(Y)[2] * dim(Y)[1])
)
nox_df <- cbind(nox_df,c(Y))
names(nox_df)[4]<-"EM_NO2"
nox_df$Longitude <- as.numeric(nox_df$Longitude)
nox_df$Latitude <- as.numeric(nox_df$Latitude)
nox_df$time <- as.Date(as.numeric(nox_df$time))
nox_df <- nox_df[order(nox_df$time,nox_df$Latitude,nox_df$Longitude),]
nox_sp <- unique(nox_df[,c("Longitude","Latitude")])
coordinates(nox_sp) <- c("Longitude","Latitude")
gridded(nox_sp)<-TRUE
nox_st <- STFDF(nox_sp, unique(nox_df$time),nox_df)
spplot(nox_st[,1,"EM_NO2"])
rm(nox_df,nox_sp)

crs_wgs84 <- CRS(SRS_string = "EPSG:4326") #reference system for WGS 84
slot(cams_st@sp, "proj4string") <- crs_wgs84
slot(we_st_land@sp, "proj4string") <- crs_wgs84
slot(we_st_sl@sp, "proj4string") <- crs_wgs84
slot(nox_st@sp, "proj4string") <- crs_wgs84

cams_era5land <- over(cams_st,we_st_land)
cams_era5sl <- over(cams_st,we_st_sl)
cams_em <- over(cams_st,nox_st)

#altitude
A5b_la_dem <- raster("data/input/LA_DEM_v100_raster_extended.tif")
A5b_la_dem <- rast(A5b_la_dem)
gridded(cams_st@sp)<-T
CAMS_res <- cams_st@sp@grid@cellsize[1]
upscal <- unique(round(CAMS_res / res(A5b_la_dem),0))
A5b_la_dem <- aggregate(A5b_la_dem,upscal,"mean")
# plot(A5b_la_dem)
# gridded(cams_st@sp)<-F
# cams_dem <- extract(raster(A5b_la_dem), cams_st@sp)
# cams_st@data <- cbind(cams_st@data,cams_dem2)
# spplot(cams_st[,1,"cams_dem2"])

cams_era5land <- cbind(cams_st@data,cams_era5land)
cams_era5sl <- cbind(cams_st@data,cams_era5sl)
cams_em <- cbind(cams_st@data,cams_em)
cams_dem <- cbind(cams_st@data,cams_dem)

names(cams_era5sl)<-paste0("sl_",names(cams_era5sl))
cams_we <- cbind(cams_era5land,cams_era5sl[,-c(1:7)])
for (i in 8:15) {
  var <- names(cams_we)[i]
  n_idx <- grep(paste0("sl_",var),names(cams_we))
  cams_we[is.na(cams_we[,i]),i] <- cams_we[is.na(cams_we[,i]),n_idx]
}

cams_we <- cams_we[,c(1:4,8:13,15:16)]
cams_we_em <- cbind(cams_we,cams_em[,"EM_NO2"])
names(cams_we_em)[ncol(cams_we_em)]<-"EM_NO2"
cams_points <- length(cams_st@sp)
cams_we_em$lag_cams_no2 <- c(rep(NA,nrow(unique(cams_we_em[,c("Longitude",
                                                              "Latitude")]))),cams_we_em$NO2[-c((nrow(cams_we_em)-(cams_points-1)):nrow(cams_we_em))])

cams_we_em <- cbind(cams_we_em,cams_dem[,"cams_dem"])
names(cams_we_em)[ncol(cams_we_em)] <- "DEM"
cams_we_em$EM_NO2 <- cams_we_em$EM_NO2*(10^6)*(86400)
summary(cams_we_em)
summary(cams_we_em[!is.na(cams_we_em$EM_NO2),])

#check
cams_sp <- unique(cams_we_em[,c("Longitude","Latitude")])
coordinates(cams_sp)<-c("Longitude","Latitude")
gridded(cams_sp)<-TRUE
cams_spdf <- SpatialPixelsDataFrame(cams_sp,data=cams_we_em[cams_we_em$time==sample(cams_we_em$time,1),])
for (i in 1:ncol(cams_we_em)) {
  print(i)
  if(i %in% c(3)){next}
  print(spplot(cams_we_em,names(cams_spdf@data)[i],main=names(cams_spdf@data)[i]))
}
cams_we_em <- cams_we_em[!is.na(cams_we_em$lai_hv),]
saveRDS(cams_we_em,file = "data/cams_we_em_2023augmented.rds")

# # ARX1 - CAMS ####
# cams_we_em$Longitude <- as.numeric(cams_we_em$Longitude)
# cams_we_em$Latitude <- as.numeric(cams_we_em$Latitude)
# lm0 <- lm(NO2 ~ .-Longitude - Latitude -time,data=cams_we_em)
# summary(lm0)
# lm1 <- lm(NO2 ~ .-Longitude - Latitude -time -lag_cams_no2 ,data=cams_we_em)
# summary(lm1)
# lm2 <- lm(NO2 ~ .-time ,data=cams_we_em)
# summary(lm2)
# 
# # without CAMS -> R2 di 0.86 a 0.4
# # remember we are using CAMS as target variable
# #
# cams_we_em_ext <- cams_we_em
# cams_we_em_ext$predict_lm <- lm0$fitted.values
# cams_we_em_ext$residuals_ls <- cams_we_em_ext$NO2 - cams_we_em_ext$predict_lm
# cams_we_em_ext$predict_lm1 <- lm1$fitted.values
# cams_we_em_ext$residuals_ls1 <- cams_we_em_ext$NO2 - cams_we_em_ext$predict_lm1
# cams_we_em_ext$predict_lm2 <- lm2$fitted.values
# cams_we_em_ext$residuals_ls2 <- cams_we_em_ext$NO2 - cams_we_em_ext$predict_lm2
# names(cams_we_em_ext)[grep("residuals",names(cams_we_em_ext))] <- paste0("res_",c("ARX1","LM","ARX1S"))
# AQ_CAMS_v100_ST@data <- cams_we_em_ext
# stplot(AQ_CAMS_v100_ST[,1:2,c("res_LM")],main="LM")
# stplot(AQ_CAMS_v100_ST[,1:2,c("res_ARX1")],main="ARX1")
# stplot(AQ_CAMS_v100_ST[,1:2,c("res_ARX1S")],main="ARX1S")
# 
# spplot(AQ_CAMS_v100_ST[,1],c("res_LM","res_ARX1","res_ARX1S"),layout=c(3,1))
# spplot(AQ_CAMS_v100_ST[,1],c("res_ARX1","res_ARX1S"),layout=c(2,1))
# 
# spplot(AQ_CAMS_v100_ST[,5],c("res_LM","res_ARX1","res_ARX1S"),layout=c(3,1))
# spplot(AQ_CAMS_v100_ST[,5],c("res_ARX1","res_ARX1S"),layout=c(2,1))
# 
# par(mfrow=c(2,1))
# stplot(AQ_CAMS_v100_ST[,4:7,c("res_LM")],main="LM")
# stplot(AQ_CAMS_v100_ST[,4:7,c("res_ARX1")],main="ARX1")
# par(mfrow=c(1,1))
# 
# saveRDS(cams_we_em_ext,file = "data/cams_we_em_2023.rds")
# 
# # preparing the prediction grid ####
cams_we_em_ext <- readRDS("data/cams_we_em_2023.rds")
crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
pred_grid_res <- 0.05
pred_grid_sp <- expand.grid(seq(min(cams_we_em_ext$Longitude),
                                max(cams_we_em_ext$Longitude),
                                pred_grid_res),
                            seq(min(cams_we_em_ext$Latitude),
                                max(cams_we_em_ext$Latitude),
                                pred_grid_res))
names(pred_grid_sp)<-c("Longitude","Latitude")
coordinates(pred_grid_sp) <- c("Longitude","Latitude")
pred_grid_st <- STFDF(sp=pred_grid_sp,
                      time=unique(as.Date(cams_we_em_ext$time)),
                      data=as.data.frame(1:(length(pred_grid_sp)*length(unique(cams_we_em_ext$time)))))

slot(pred_grid_st@sp, "proj4string") <- crs_wgs84
pred_grid_df <- over(pred_grid_st,AQ_CAMS_v100_ST)

A5a_la_clc <- raster("data/input/LA_CLC_v100_raster.tif")
A5a_la_clc <- rast(A5a_la_clc)
BAUs_res <- c(BAUs_res,BAUs_res)
upscal <- unique(round(BAUs_res / res(A5a_la_clc),0))
A5a_la_clc <- aggregate(A5a_la_clc,upscal,"modal")
pred_grid_df$CLC <- extract(raster(A5a_la_clc), pred_grid_st@sp)

A5b_la_dem <- raster("data/input/LA_DEM_v100_raster.tif")
A5b_la_dem <- rast(A5b_la_dem)
upscal <- unique(round(BAUs_res / res(A5b_la_dem),0))
A5b_la_dem <- aggregate(A5b_la_dem,upscal,"mean")
pred_grid_df$DEM <- extract(raster(A5b_la_dem), pred_grid_st@sp)

# plot(A5a_la_clc,type = "classes")
# plot(A5b_la_dem)
names(pred_grid_df) <- gsub("Longitude","Lon_CAMS",names(pred_grid_df))
names(pred_grid_df) <- gsub("Latitude","Lat_CAMS",names(pred_grid_df))
pred_grid_df <- cbind(pred_grid_sp@coords,pred_grid_df)
saveRDS(pred_grid_df,file = "data/pred_grid_df.rds")
# pred_grid_df <- readRDS("data/pred_grid_df.rds")


# AQ-EEA ####
load("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-EEA/v.1.0.3/data/Zenodo/GRINS_AQCLIM_points_Italy.rda")
aq_eea <- GRINS_AQCLIM_points_Italy[format(GRINS_AQCLIM_points_Italy$time,"%Y") == "2023",]
rm(GRINS_AQCLIM_points_Italy)
names(aq_eea)
aq_eea <- aq_eea[order(aq_eea$time,aq_eea$AirQualityStation),]
length(unique(aq_eea$AirQualityStation))

load("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-EEA/v.1.0.5/data/Zenodo/Station_registry_information.rda")
names(Station_registry_information)
length(unique(Station_registry_information$AirQualityStation))
setdiff(unique(aq_eea$AirQualityStation),unique(Station_registry_information$AirQualityStation))
meta_aq_eea <- aq_eea[aq_eea$time == aq_eea$time[1],]
meta_aq_eea <- merge(meta_aq_eea,Station_registry_information)
sp_aq_eea <- meta_aq_eea[,c("AirQualityStation","Longitude","Latitude")]
names(sp_aq_eea) <- gsub("AirQualityStation","Station_ID",names(sp_aq_eea))
coordinates(sp_aq_eea) <- c("Longitude","Latitude")
st_aq_eea <- STFDF(sp_aq_eea,unique(aq_eea$time),aq_eea)

pred_grid_df <- readRDS("data/pred_grid_df.rds")
pred_grid_df <- pred_grid_df[order(pred_grid_df$time,pred_grid_df$Latitude,pred_grid_df$Longitude),]
meta_pred <- pred_grid_df[,c("Longitude","Latitude")]
meta_pred <- unique(meta_pred)
coordinates(meta_pred) <- c("Longitude","Latitude")
gridded(meta_pred)<-TRUE
st_pred <- STFDF(meta_pred,as.Date(unique(pred_grid_df$time)),pred_grid_df)

aq_eea_ext <- cbind(st_aq_eea@data,over(st_aq_eea,st_pred))

names(aq_eea_ext)
aq_eea_df <- aq_eea_ext[,c("AirQualityStation",
              "time",
              "AQ_mean_NO2",names(pred_grid_df)[-c(1:5)])]
names(aq_eea_df)[c(3,4)]<-c("EEA_NO2","CAMS_NO2")
saveRDS(aq_eea_df,file = "data/aq_eea_df.rds")
