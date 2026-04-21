# weather from matrix to ST ####
library(spacetime)
library(abind)
# ear5land ####
lf <- list.files("data/ext_22_24/we_2024","rds")
lf <- lf[-which(lf=="blh.rds")]
ndata <- length(lf)
for (i in 1:ndata) {
  y <- readRDS(paste0("data/ext_22_24/we_2024/", lf[i]))
  if (i == 1) {
    Y <- y
  } else{
    Y <- abind(Y, y, along = 4)
  }
}
varnames <- gsub(".rds","",lf)
dimnames(Y)[[4]] <- varnames
rm(y)
WE_C3S_v100_df_ERA5Land <- data.frame(
  Longitude = rep(rep(dimnames(Y)[[2]], each = dim(Y)[1]), dim(Y)[3]),
  Latitude = rep(dimnames(Y)[[1]], dim(Y)[2] * dim(Y)[3]),
  time = rep(as.Date(as.numeric(dimnames(
    Y
  )[[3]])), each = dim(Y)[2] * dim(Y)[1])
)
WE_C3S_v100_df_ERA5Land <- cbind(WE_C3S_v100_df_ERA5Land,
                                 matrix(c(Y), ncol = dim(Y)[4]))
names(WE_C3S_v100_df_ERA5Land)[-c(1:3)] <- varnames
centre <-
  data.frame(Longitude = as.numeric(rep(rep(
    dimnames(Y)[[2]], each = dim(Y)[1]
  ))),
  Latitude = as.numeric(rep(dimnames(Y)[[1]], dim(Y)[2])))
coordinates(centre) <- c("Longitude", "Latitude")
gridded(centre) <- TRUE
colnames(centre@coords) <- c("coords.x1", "coords.x2")
WE_C3S_v100_ST_ERA5Land <- STFDF(
  sp = centre,
  time = unique(WE_C3S_v100_df_ERA5Land$time),
  data = WE_C3S_v100_df_ERA5Land
)
# png(file="AQ-FRK/v.3.0.0/A_input/plot/input/WE_C3S_ERA5Land_day1e2.png")
# stplot(WE_C3S_v100_ST_ERA5Land[, 1:2, "ssr" ])
# dev.off()

save(WE_C3S_v100_ST_ERA5Land,
     file = "AQ-FRK/v.3.0.0/A_input/data/input/WE_C3S_v100_ST_ERA5Land.rda")
rm(list=setdiff(ls(),c(start_ls,"start_ls")))
gc()

# A3b_we_era5sl ####
lf <- list.files("WE-C3S/v.1.0.0/data/ERA5SL/daily","Rdata")
lf <- lf[-8]

ndata <- length(lf)
for (i in 1:ndata) {
  if(lf[i]=="daily_tp.Rdata"){
    load(paste0("WE-C3S/v.1.0.1/data/ERA5SL/daily/", lf[i]))
  }else{
    load(paste0("WE-C3S/v.1.0.0/data/ERA5SL/daily/", lf[i]))
  }
  if (i == 1) {
    Y <- y
  } else{
    Y <- abind(Y, y, along = 4)
  }
}
varnames <- unlist(lapply(lf, function(x)
  substr(x, 7, nchar(x) - 6)))
dimnames(Y)[[4]] <- varnames
rm(y)
WE_C3S_v100_df_ERA5SL <- data.frame(
  Longitude = rep(rep(dimnames(Y)[[2]], each = dim(Y)[1]), dim(Y)[3]),
  Latitude = rep(dimnames(Y)[[1]], dim(Y)[2] * dim(Y)[3]),
  time = rep(as.Date(as.numeric(dimnames(
    Y
  )[[3]])), each = dim(Y)[2] * dim(Y)[1])
)
WE_C3S_v100_df_ERA5SL <- cbind(WE_C3S_v100_df_ERA5SL,
                               matrix(c(Y), ncol = dim(Y)[4]))
names(WE_C3S_v100_df_ERA5SL)[-c(1:3)] <- varnames
centre <-
  data.frame(Longitude = as.numeric(rep(rep(
    dimnames(Y)[[2]], each = dim(Y)[1]
  ))),
  Latitude = as.numeric(rep(dimnames(Y)[[1]], dim(Y)[2])))
coordinates(centre) <- c("Longitude", "Latitude")
gridded(centre) <- TRUE
colnames(centre@coords) <- c("coords.x1", "coords.x2")
WE_C3S_v100_ST_ERA5SL <- STFDF(
  sp = centre,
  time = unique(WE_C3S_v100_df_ERA5SL$time),
  data = WE_C3S_v100_df_ERA5SL
)
png(file="AQ-FRK/v.3.0.0/A_input/plot/input/WE_C3S_ERA5SL_day1e2.png")
stplot(WE_C3S_v100_ST_ERA5SL[, 1:2, "blh"])
dev.off()
save(WE_C3S_v100_ST_ERA5SL,
     file = "AQ-FRK/v.3.0.0/A_input/data/input/WE_C3S_v100_ST_ERA5SL.rda")
rm(list=setdiff(ls(),c(start_ls,"start_ls")))
gc()

# emissions
#

# from preparing datasets
load("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-CAMS/v.1.0.0/data/output/ST/AQ_CAMS_v100_ST_no2.rda")
AQ_CAMS_v100_ST_0 <- AQ_CAMS_v100_ST[,which(index(AQ_CAMS_v100_ST@time)== t[1]-1)]
AQ_CAMS_v100_ST <- AQ_CAMS_v100_ST[,which(index(AQ_CAMS_v100_ST@time)%in% t)]
load("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-FRK/v.3.0.0/A_input/data/input/WE_C3S_v100_ST_ERA5SL.rda")
WE_C3S_v100_ST_ERA5SL <- WE_C3S_v100_ST_ERA5SL[,which(index(WE_C3S_v100_ST_ERA5SL@time)%in% t)]

load("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-FRK/v.3.0.0/A_input/data/input/WE_C3S_v100_ST_ERA5Land.rda")
WE_C3S_v100_ST_ERA5Land <- WE_C3S_v100_ST_ERA5Land[,which(index(WE_C3S_v100_ST_ERA5Land@time)%in% t)]
EM_CAMS_v300_ST <- readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/EM-CAMS/v.3.0.0/data.nosync/ST/EM_NO2_CAMS_v300_ST_2023.rds")
