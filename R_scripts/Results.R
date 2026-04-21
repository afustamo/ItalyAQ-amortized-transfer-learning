library(sf)
library(ggplot2)

cv_f <- list.files("data/output/CV",pattern = ".rds")
i <- cv_f[1]
Y <- readRDS(file.path("data/output/CV",i))
for (i in cv_f[-1]) {
  y <- readRDS(file.path("data/output/CV",i))
  Y <- rbind(Y,y)
}

write.csv(Y,"data/output/CV_results.csv")

eea_df_2023 <- readRDS("data/eea_df_2023.rds")

cv_df <- merge(eea_df_2023,Y)
length(unique(cv_df$time))

summary(cv_df)

# printing overall results ####
rmse<-c()
for (i in grep("_fit",names(cv_df))){
  print(names(cv_df)[i])
  rmse <- c(rmse,sqrt(mean((cv_df$EEA_NO2 - cv_df[,i])^2)))
}
rmse
# rmse
# [1] 9.967256 9.940797 9.684350 9.736979 9.661096
# 0.3 / 9.967256
# [1] 0.03009855

picp <- c()
for (i in grep("_lwr",names(cv_df))){
  print(names(cv_df)[i])
  cpi <- cv_df$EEA_NO2 >= cv_df[,i] & cv_df$EEA_NO2 <= cv_df[,i+2]
  picp <- c(picp,sum(cpi)/length(cpi))
}
picp

mpiw <- c()
for (i in grep("_lwr",names(cv_df))){
  print(names(cv_df)[i])
  piw <- cv_df[,i+2] - cv_df[,i]
  mpiw <- c(mpiw,mean(piw))
}
mpiw

# plotting time series overall error by categories ####

cv_df$e_lm <- cv_df$EEA_NO2 - cv_df$lm_fit
cv_df$e_stat <- cv_df$EEA_NO2 - cv_df$LK_stat_fit
cv_df$e_stun <- cv_df$EEA_NO2 - cv_df$LK_stun_adj_fit

library(ggplot2)

start <- as.Date("2023-02-01")
end <- as.Date("2023-02-28")
cv_df_sub <- cv_df[cv_df$AirQualityStation == "STA.IT0502A",]
ggplot(cv_df_sub,aes(x=time))+
  geom_hline(yintercept = 0,linetype=2)+
  geom_line(aes(y=e_lm),col="green",alpha=.5)+
  geom_line(aes(y=e_stat),col="blue",,alpha=.5)+
  geom_line(aes(y=e_stun),col="orange",alpha=.5)+
  coord_cartesian(xlim = c(start,end))+
  ggtitle(label=unique(cv_df_sub$AirQualityStation))

library(readr)
Station_registry_information <- read_csv("data/Station_registry_information.CSV")
cv_df <- merge(cv_df,Station_registry_information[,c("AirQualityStation","AirQualityStationType","AirQualityStationArea")])
# cv_df$AirQualityStationArea
# unique(cv_df$AirQualityStationType)
# unique(cv_df$AirQualityStationArea)
meta_cv <- unique(cv_df[,c("AirQualityStation","AirQualityStationType","AirQualityStationArea")])
table(meta_cv[,c("AirQualityStationType","AirQualityStationArea")])
cv_df$AirQualityStationArea <- gsub("-.*","",cv_df$AirQualityStationArea)
# unique(cv_df$AirQualityStationArea)
ggplot(cv_df)+
  geom_boxplot(aes(y=e_lm,col=AirQualityStationType))+
  facet_wrap(~AirQualityStationArea)
  # geom_boxplot(aes(y=e_stat))

ggplot(cv_df)+
  geom_boxplot(aes(y=e_stat,col=AirQualityStationType))+
  facet_wrap(~AirQualityStationArea)


ggplot(cv_df)+
  geom_boxplot(aes(y=e_stun,col=AirQualityStationType))+
  facet_wrap(~AirQualityStationArea)

t.test(cv_df$e_stat,
       cv_df$e_stun,
       paired = T)

# Map of errors by province ####
# cv_df$diff <- abs(cv_df$e_stat) - abs(cv_df$e_stun)
cv_df$e_stat2 <- cv_df$e_stat^2
cv_df$e_stun2 <- cv_df$e_stun^2

# summary(cv_df$diff)
# hist(cv_df$diff[cv_df$time == as.Date("2023/02/13")])
# # cv_df$diff <- cv_df$diff / cv_df$e_stat
# cv_df$diff2 <- abs(cv_df$e_stun) - abs(cv_df$e_lm)
# 
# a <- cv_df[order(cv_df$diff),]
# b <- cv_df[order(cv_df$diff2),]
# summary(a)


library(dplyr)

cv_df_feb <- cv_df[months(cv_df$time)=="February",]
cv_df_feb <- cv_df

aa <- aggregate(cv_df_feb$e_stat2,list(cv_df_feb$AirQualityStation),"mean")
bb <- aggregate(cv_df_feb$e_stun2,list(cv_df_feb$AirQualityStation),"mean")

aa$x <- sqrt(aa$x)
bb$x <- sqrt(bb$x)
names(aa)<-c("AirQualityStation","RMSE")
names(bb)<-c("AirQualityStation","RMSE")
aa <- merge(aa,Station_registry_information)
bb <- merge(bb,Station_registry_information)
aaa <- aggregate(aa$RMSE,list(aa$COD_PROV),"mean")
bbb <- aggregate(bb$RMSE,list(bb$COD_PROV),"mean")
names(aaa)<-c("COD_PROV","RMSE")
names(bbb)<-c("COD_PROV","RMSE")
ccc <- cbind(aaa,bbb[,2])
names(ccc)[3]<-"RMSE_stun"
wld_shp <- st_read("data/shapefile/CNTR_RG_60M_2024_4326")

xlim <- c(6.0,18.9)
ylim <- c(35.0,47.9)

mun <- geotools::prov_bounds
mun <- st_transform(mun,st_crs(wld_shp))
aaaa <- merge(mun,ccc,all.x=T)
summary(aaaa$RMSE)
library(scales)

# ggplot()+
#   geom_sf(data = aaaa,aes(fill=diff))+
#   scale_fill_gradient2(
#     low = muted("red"),
#     mid = "white",
#     high = muted("darkgreen"),
#     midpoint = 0,
#     limits = c(-6, 6),
#     trans = pseudo_log_trans(sigma = 0.02),
#     name = "Diff"
#   )+ geom_sf(
#     data = wld_shp,
#     fill = NA,
#     linewidth = .4,
#     col = "black"
#   ) 

ggplot()+
  geom_sf(
    data = wld_shp,
    fill = NA,
    linewidth = .4,
    col = "black"
  ) 
  
# ggplot()+
#   geom_sf(data = aaaa,aes(fill=diff))

st_crs(aaaa)
st_crs(wld_shp)

png("plot/RMSEbyprovince.png",width = 300, height = 300)
ggplot()+
  geom_sf(
    data = wld_shp[!wld_shp$CNTR_ID=="IT",],
    fill = NA,
    linewidth = .4,
    col = "black"
  ) +
  geom_sf(data = aaaa,aes(fill=RMSE_stun - RMSE))+
  coord_sf(
    xlim = xlim,
    ylim = ylim,
    expand = FALSE
    # label_axes = "E--N"
  )+
  scale_fill_gradient2(
        low = muted("darkgreen"),
        mid = "white",
        high = muted("red"),
        midpoint = 0,
        limits = c(-3, 3),
        trans = pseudo_log_trans(sigma = 0.02),
        name = "RMSE_stun -
RMSE_stat",
        breaks = c(-3,-1,-0.25,0,0.25,1,3)
      )
dev.off()
  
  # scale_fill_continuous(transform="pseudo_log",
  #                       type="viridis")+
  # scale_fill_ma

# matplot()

cv_df_time <- aggregate(cv_df$e_stat - cv_df$e_stun,list(cv_df$time),mean)
names(cv_df_time)<-c("time","err")

ggplot(cv_df_time)+
  geom_point(aes(y=err,x=time))+
  geom_hline(yintercept = 0,linetype=1)+
  geom_smooth(
    aes(x = time, y = err),
    method = "gam",
    formula = y ~ s(x, bs = "tp",k=10),
    se = TRUE
  )

# cfr kappa
library(ggplot2)
LK_stat <- readRDS("data/output/LK_stat_14.rds")
LK_stun_adj <- readRDS("data/output/LK_stun_adj_14.rds")
a<-LK_stun_adj[["LKinfo"]][["a.wghtObject"]][["z"]][,,5] + 
  2*LK_stun_adj[["LKinfo"]][["a.wghtObject"]][["z"]][,,2] + 
  2*LK_stun_adj[["LKinfo"]][["a.wghtObject"]][["z"]][,,4]
ax<-LK_stun_adj[["LKinfo"]][["a.wghtObject"]][["x"]]
ay<-LK_stun_adj[["LKinfo"]][["a.wghtObject"]][["y"]]

image(a)
# a <- LK_stun_adj[["MLE"]][["LKinfo"]][["a.wghtObject"]][["z"]][,,5] - 4
# ax <- LK_stun_adj[["MLE"]][["LKinfo"]][["a.wghtObject"]][["x"]]
# ay <- LK_stun_adj[["MLE"]][["LKinfo"]][["a.wghtObject"]][["y"]]

# head(expand.grid(ax,ay))
b <- LK_stat[["MLE"]][["a.wght.MLE"]] - 4
bb <- matrix(b,dim(a),dim(a))
df <- expand.grid(ax,ay)
df$y <- c(a)
df$y2 <- c(b)
ggplot(df)+
  geom_tile(aes(Var1,Var2,fill=y-y2))

# guardando al kappa del ns mi sembra che sia molto elevato dove non dovrebbe...
# anche se effettivamente nel mare è molto scuro...
# bisognerebbe metterci l'Italia sopra...
library(sf)
library(sp)
coordinates(df)<-c("Var1","Var2")
gridded(df)<-TRUE
italy_plot <- spplot(Italy)
a_plot <- spplot(df,"y")
a_plot + italy_plot
load(
  "/Users/alessandrofusta/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-ALL/v.1.0.0/data/Italy_shp.Rdata"
)
Italy <-st_as_sf(Italy)
st_crs(Italy)
# Italy <- as_Spatial(Italy)
# plot(Italy,add=T)

df_sf <- st_as_sf(df)
st_crs(df_sf) <- st_crs(4326)
ggplot()+
  geom_sf(data = df_sf,aes(col=y))+
  geom_sf(data = Italy,fill=NA,col="black",linewidth=1)


