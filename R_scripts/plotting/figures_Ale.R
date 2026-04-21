library(sp)
library(sf)
library(readr)
library(patchwork)

# cams_we_em_2023 <- readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/Italy_AQ_AmortisedLatticeKrig/data/cams_we_em_2023.rds")
aq_eea_df <- readRDS(
  "~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/Italy_AQ_AmortisedLatticeKrig/data/aq_eea_df.rds"
)
pred_grid_df <- readRDS(
  "~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/Italy_AQ_AmortisedLatticeKrig/data/pred_grid_df.rds"
)

t_sel <- unique(aq_eea_df$time)[103]
sub <- pred_grid_df[pred_grid_df$time == t_sel, ]
sub_eea <- aq_eea_df[aq_eea_df$time == t_sel, ]
Station_registry_information <- read_csv("data/Station_registry_information.CSV")
sub_eea <- merge(sub_eea, Station_registry_information[, c("AirQualityStation", "Longitude", "Latitude")], all.x =
                   T)
pred_sp <- unique(sub[, c("Longitude", "Latitude")])

coordinates(pred_sp) <- c("Longitude", "Latitude")
pred_sf <- st_as_sf(pred_sp)
st_crs(pred_sf) <- st_crs(4326)

load(
  "/Users/alessandrofusta/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2/AQ-ALL/v.1.0.0/data/Italy_shp.Rdata"
)
Italy <- st_as_sf(Italy)
Italy <- st_transform(Italy, st_crs(pred_sf))

# add world shp
wld_shp <- st_read("data/shapefile/CNTR_RG_60M_2024_4326")
# ggplot()+
#   geom_sf(data=wld_shp)

library(ggplot2)
xlim <- range(sub$Longitude, na.rm = TRUE)
ylim <- range(sub$Latitude, na.rm = TRUE)

col_C_lims <- c(min(sub_eea$EEA_NO2[!is.na(sub_eea$EEA_NO2)]),
                  max(sub_eea$EEA_NO2[!is.na(sub_eea$EEA_NO2)]))

breaks_c_col <- c(0,1,2,5,10,15,20,25,35,50)

p4 <- ggplot() +
  geom_point(data = sub_eea[!is.na(sub_eea$EEA_NO2), ], aes(Longitude, Latitude, col = EEA_NO2)) +
  scale_color_viridis_c(option = "C", name = expression("µg/m"^3), limits=col_C_lims, 
                        transform="pseudo_log", breaks=breaks_c_col) +
  geom_sf(
    data = wld_shp,
    fill = NA,
    linewidth = .4,
    col = "black"
  ) +
  coord_sf(
    xlim = xlim,
    ylim = ylim,
    expand = FALSE,
    label_axes = "E--N"
  ) +
  theme(axis.title.x.bottom = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.title.position = "bottom",
        legend.title.align = 0.5)+
  scale_x_continuous(breaks = seq(xlim[1]+2,xlim[2], by = 3)) +
  scale_y_continuous(breaks = seq(ylim[1],ylim[2], by = 3)) + guides(
    color = guide_colorbar(
      title.position = "bottom",
      barwidth = unit(10, "cm")  # allunga la barra
      # barheight = unit(0.5, "cm")
    )
  )
  # scale_y_continuous(breaks = c(6, 10, 14)) 


p1 <- ggplot() +
  geom_tile(data = sub, aes(Longitude, Latitude, fill = NO2)) +
  scale_fill_viridis_c(option = "C", name = expression("µg/m"^3), limits=col_C_lims, 
                       transform="pseudo_log", breaks=breaks_c_col) +
  geom_sf(
    data = wld_shp,
    fill = NA,
    linewidth = .4,
    col = "black"
  ) +
  coord_sf(
    xlim = xlim,
    ylim = ylim,
    expand = FALSE,
    label_axes = "E--N"
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        legend.position = "bottom",
        legend.title.position = "bottom",
        legend.title.align = 0.5)+
  scale_x_continuous(breaks = seq(xlim[1]+2,xlim[2], by = 3)) +
  scale_y_continuous(breaks = seq(ylim[1],ylim[2], by = 3)) + guides(
    fill = guide_colorbar(
      title.position = "bottom",
      barwidth = unit(10, "cm")  # allunga la barra
      # barheight = unit(0.5, "cm")
    )
  )

p2 <- ggplot() +
  geom_tile(data = sub, aes(Longitude, Latitude, fill = t2m)) +
  scale_fill_distiller(
    palette = "Spectral",
    name = expression(degree * C)
  ) +
  geom_sf(
    data = wld_shp,
    fill = NA,
    linewidth = .4,
    col = "black"
  )+
  coord_sf(
    xlim = xlim,
    ylim = ylim,
    expand = FALSE,
    label_axes = "E--N"
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        legend.position = "bottom",
        legend.title.position = "bottom",
        legend.title.align = 0.5)+
  scale_x_continuous(breaks = seq(xlim[1]+2,xlim[2], by = 3)) +
  scale_y_continuous(breaks = seq(ylim[1],ylim[2], by = 3)) 
# labs(title = "Temperatura")

p3 <- ggplot() +
  geom_tile(data = sub, aes(Longitude, Latitude, fill = EM_NO2)) +
  scale_fill_viridis_c(option = "E", name = expression("mg/m"^2)) +
  geom_sf(
    data = wld_shp,
    fill = NA,
    linewidth = .4,
    col = "black"
  )+
  coord_sf(
    xlim = xlim,
    ylim = ylim,
    expand = FALSE,
    label_axes = "E--N"
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        legend.position = "bottom",
        legend.title.position = "bottom",
        legend.title.align = 0.5)+
  scale_x_continuous(breaks = seq(xlim[1]+2,xlim[2], by = 3)) +
  scale_y_continuous(breaks = seq(ylim[1],ylim[2], by = 3)) 
# labs(title = expression("Emission of "*NO[2]))


AQ_plot <- p4 | p1
AQ_plot <- AQ_plot + plot_layout(guides = "collect")
AQ_plot <- AQ_plot &
  theme(legend.position = "bottom",
        legend.box.just = "center",
        legend.margin = margin(r = 100,l=-5))

p2 <- p2 &
  theme(legend.margin = margin(r = 65,l=-55))

# 
# AQ_plot
                                 
(AQ_plot | p2 | p3) 
# &
#   theme(legend.position = "bottom",legend.box.just = "center")
#     &
#   theme(
#     # legend.position = "bottom",
#     legend.box.spacing = unit(0.2, "cm"),
#     # legend.margin = margin(t = 2,r=20,l=20),
#     legend.title.position = "bottom",
#     legend.title.align = 0.5
#   )

ggsave("plot/data_input.pdf", width = 10, height = 4.3)

# 
# theme_nox <- theme(
#   axis.title.x = element_blank(),
#   axis.text.x  = element_blank(),
#   axis.ticks.x = element_blank()
# )
# 
# theme_noy <- theme(
#   axis.title.y = element_blank(),
#   axis.text.y  = element_blank(),
#   axis.ticks.y = element_blank()
# )
# 
# # p1 | p2
# # p3 | p4
# # label_axes = list(bottom = "E", top="E", left="N", right="N")
# 
# p4 <- p4 +
#   coord_sf(
#     xlim = xlim,
#     ylim = ylim,
#     expand = FALSE,
#     label_axes = "E--N"
#   )
# # theme(
# #   axis.title.x.bottom = element_blank(),
# #   axis.text.x.bottom  = element_blank(),
# #   axis.ticks.x.bottom = element_blank()
# # )
# p1 <- p1 +
#   coord_sf(xlim = xlim,
#            ylim = ylim,
#            expand = FALSE) +
#   theme_noy +
#   scale_x_continuous(position = "top")
# p2 <- p2 +
#   coord_sf(xlim = xlim,
#            ylim = ylim,
#            expand = FALSE) +
#   theme_noy +
#   scale_x_continuous(position = "top")
# p3 <- p3 +
#   coord_sf(xlim = xlim,
#            ylim = ylim,
#            expand = FALSE) +
#   theme_noy +
#   scale_x_continuous(position = "top")
# 
# # library(patchwork)
# #
# # (p1 | p2) /
# #   (p3 | p4)
# 
# 
# 
# # pred_spp <- SpatialPointsDataFrame(pred_sp,pred_grid_df[pred_grid_df$time==t_sel,])
# # gridded(pred_spp) <- TRUE
# # spplot(pred_spp,c("NO2"))
# # spplot(pred_spp,c("t2m"))
# # spplot(pred_spp,c("EM_NO2"))
# # spplot(pred_spp,c("DEM"))

# # μg/m3
# ggplot() +
#   geom_tile(data = sub, aes(Longitude, Latitude, fill = NO2)) +
#   scale_fill_continuous(type = "viridis") +
#   geom_sf(
#     data = Italy,
#     fill = NA,
#     linewidth = .5,
#     col = "black"
#   ) +
#   
#   # °C
#   ggplot() +
#   geom_tile(data = sub, aes(Longitude, Latitude, fill = t2m)) +
#   scale_fill_continuous(type = "viridis") +
#   geom_sf(
#     data = Italy,
#     fill = NA,
#     linewidth = .5,
#     col = "black"
#   )
# 
# # mg/m2
# ggplot() +
#   geom_tile(data = sub, aes(Longitude, Latitude, fill = EM_NO2)) +
#   scale_fill_continuous(type = "viridis") +
#   geom_sf(
#     data = Italy,
#     fill = NA,
#     linewidth = .5,
#     col = "black"
#   )
# 
# # m
# ggplot() +
#   geom_tile(data = sub, aes(Longitude, Latitude, fill = DEM)) +
#   scale_fill_continuous(type = "viridis") +
#   geom_sf(
#     data = Italy,
#     fill = NA,
#     linewidth = .5,
#     col = "black"
#   )

