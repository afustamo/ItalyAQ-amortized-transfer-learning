# cv_f <- list.files("data/output/CV")
# i <- cv_f[1]
# Y <- readRDS(file.path("data/output/CV",i))
# for (i in cv_f[-1]) {
#   y <- readRDS(file.path("data/output/CV",i))
#   Y <- rbind(Y,y)
# }
# import hacked LK files for additional kappa2 estimation
setwd("C:/Users/anton/Desktop/Research/Italy_AQ_paper")

# Y <- read.csv("Italy_AQ_AmortisedLatticeKrig/results/full_crossval_results.csv")
Y <- read.csv("Italy_AQ_AmortisedLatticeKrig/results/10fold_results.csv")
eea_df_2023 <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/eea_df_2023.rds")

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


library(dplyr)

# helper to compute coverage given lwr/upr column names
station_picp_one <- function(df, lwr_col, upr_col) {
  df %>%
    mutate(covered = .data$EEA_NO2 >= .data[[lwr_col]] & .data$EEA_NO2 <= .data[[upr_col]]) %>%
    group_by(AirQualityStation, Longitude, Latitude) %>%
    summarise(picp = mean(covered, na.rm = TRUE), .groups = "drop")
}

# per-model station-level PICP dfs
picp_lm       <- station_picp_one(cv_df, "lm_lwr",       "lm_upr")       %>% mutate(model = "lm")
picp_stat     <- station_picp_one(cv_df, "LK_stat_lwr",  "LK_stat_upr")  %>% mutate(model = "LK_stat")
picp_stun     <- station_picp_one(cv_df, "LK_stun_lwr",  "LK_stun_upr")  %>% mutate(model = "LK_stun")
picp_stun_adj <- station_picp_one(cv_df, "LK_stun_adj_lwr","LK_stun_adj_upr") %>% mutate(model = "LK_stun_adj")

# long (one row per station per model)
station_picp_long <- bind_rows(picp_lm, picp_stat, picp_stun, picp_stun_adj)

# wide (one row per station with a column per model) — usually nice for mapping
station_picp_wide <- station_picp_long %>%
  tidyr::pivot_wider(names_from = model, values_from = picp)


station_picp_wide

colormap <- colorRampPalette(brewer.pal(11, "Spectral"))(256)
# colormap <- rev(turbo(256))

# nonstat adjusted
bubblePlot(
  x = station_picp_wide$Longitude,
  y = station_picp_wide$Latitude,
  z = station_picp_wide$LK_stun_adj,
  main = "Station PICP (LK_stun_adj)",
  xlab = "Longitude",
  ylab = "Latitude",
  col = colormap
)
world(add = TRUE)



# 
# 
# 
# library(ggplot2)
# 
# start <- as.Date("2023-02-01")
# end <- as.Date("2023-02-28")
# cv_df_sub <- cv_df[cv_df$AirQualityStation == "STA.IT0502A",]
# ggplot(cv_df_sub,aes(x=time))+
#   geom_hline(yintercept = 0,linetype=2)+
#   geom_line(aes(y=e_lm),col="green",alpha=.5)+
#   geom_line(aes(y=e_stat),col="blue",,alpha=.5)+
#   geom_line(aes(y=e_stun),col="orange",alpha=.5)+
#   coord_cartesian(xlim = c(start,end))+
#   ggtitle(label=unique(cv_df_sub$AirQualityStation))
# 
# library(readr)
# Station_registry_information <- load("Italy_AQ_AmortisedLatticeKrig/data/Station_registry_information.rda")
# cv_df <- merge(cv_df,Station_registry_information[,c("AirQualityStation","AirQualityStationType","AirQualityStationArea")])
# # cv_df$AirQualityStationArea
# # unique(cv_df$AirQualityStationType)
# # unique(cv_df$AirQualityStationArea)
# meta_cv <- unique(cv_df[,c("AirQualityStation","AirQualityStationType","AirQualityStationArea")])
# table(meta_cv[,c("AirQualityStationType","AirQualityStationArea")])
# cv_df$AirQualityStationArea <- gsub("-.*","",cv_df$AirQualityStationArea)
# # unique(cv_df$AirQualityStationArea)
# ggplot(cv_df)+
#   geom_boxplot(aes(y=e_lm,col=AirQualityStationType))+
#   facet_wrap(~AirQualityStationArea)
# # geom_boxplot(aes(y=e_stat))
# 
# ggplot(cv_df)+
#   geom_boxplot(aes(y=e_stat,col=AirQualityStationType))+
#   facet_wrap(~AirQualityStationArea)
# 
# 
# ggplot(cv_df)+
#   geom_boxplot(aes(y=e_stun,col=AirQualityStationType))+
#   facet_wrap(~AirQualityStationArea)
# 
# t.test(cv_df$e_stat,
#        cv_df$e_stun,
#        paired = T)
# 
# # Map of errors by province ####
# cv_df$diff <- abs(cv_df$e_stat) - abs(cv_df$e_stun)
# summary(cv_df$diff)
# # cv_df$diff <- cv_df$diff / cv_df$e_stat
# cv_df$diff2 <- abs(cv_df$e_stun) - abs(cv_df$e_lm)
# 
# a <- cv_df[order(cv_df$diff),]
# b <- cv_df[order(cv_df$diff2),]
# summary(a)
# 
# 
# library(dplyr)
# 
# cv_df_feb <- cv_df[months(cv_df$time)=="February",]
# cv_df_feb <- cv_df
# 
# aa <- aggregate(cv_df_feb$diff,list(cv_df_feb$AirQualityStation),"mean")
# names(aa)<-c("AirQualityStation","diff")
# aa <- merge(aa,Station_registry_information)
# aaa <- aggregate(aa$diff,list(aa$COD_PROV),"mean")
# names(aaa)<-c("COD_PROV","diff")
# 
# 
# mun <- geotools::prov_bounds
# aaaa <- merge(mun,aaa,all.x=T)
# summary(aaaa$diff)
# library(scales)
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
#   )
# # scale_fill_continuous(transform="pseudo_log",
# #                       type="viridis")+
# # scale_fill_ma
# 
# # matplot()
# 
# cv_df_time <- aggregate(cv_df$e_stat - cv_df$e_stun,list(cv_df$time),mean)
# names(cv_df_time)<-c("time","err")
# 
# ggplot(cv_df_time)+
#   geom_point(aes(y=err,x=time))+
#   geom_smooth(aes(y=err,x=time))
# 
# 
# 
