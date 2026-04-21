# ------------------------------------
# Libraries
# ------------------------------------
library(rhdf5)
library(LatticeKrig)
library(sf)
library(maps)
library(cmocean)
library(scico)
library(RColorBrewer)

library(zoo) 
library(mgcv) 

# ------------------------------------
# Setup and Funcs
# ------------------------------------

# import hacked LK files for additional kappa2 estimation
setwd("C:/Users/anton/Desktop/Research/Italy_AQ_paper")

set.seed(7)

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


# CHOOSING STATIONS 
stations_df <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/eea_df_2023.rds")
cams_df     <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/cams_df_2023padded.rds")
cams_df <- cams_df[cams_df$time >= 19358 & cams_df$time <= 19722, ]



# ------------------------------------
# Param Ellipses Figures
# ------------------------------------

param_path <- "Italy_AQ_AmortisedLatticeKrig/data/STUN_param_df_2023padded.h5" #STUN_param_df_2023padded.h5
param_data_name <- "arx1_surround_30rep_output" #"arx1_surround_30rep_output"

params <- h5read(param_path, param_data_name)
H5close()
# for plotting and image purposes, need to flip upside down
params <- params[, 128:1, , ]

# recover kappa
params[,,1,] <- exp(params[,,1,])
# theta is ,,2, and rho is ,,3,

chosen_day <- 103

cams_df_day <- cams_df[
  cams_df$time == (cams_df$time[1] + chosen_day - 1),
]


# Big, scary, ugly, CoPilot function to plot ellipses 
# It has a million inputs but this is important because 
# I want to make tons of plots 
# Example uses shown below

# Ellipse drawing helper 
ellipse_xy_helper <- function(x0, y0, a, b, ang, n = 60) {
  t <- seq(0, 2*pi, length.out = n)
  x <- a * cos(t)
  y <- b * sin(t)
  xr <- x * cos(ang) - y * sin(ang)
  yr <- x * sin(ang) + y * cos(ang)
  list(x = x0 + xr, y = y0 + yr)
}

# Main function
plot_param_ellipses <- function(
    # params
  kappa2, 
  rho,
  theta,
  # your grid 
  gridList, 
  # pixel spacing between ellipses 
  pixel_spacing = 6, 
  # domain of interest
  xlim = range(gridList$x),
  ylim = range(gridList$y),
  
  # optional background field to plot
  draw_field = FALSE,
  plot_field = NULL,              # matrix, e.g. CTM NO2 on same grid
  field = "field",                # label for title/colorbar
  field_col = turbo(256),
  field_zlim = NULL,
  
  # If not NULL: exact zeros in plot_field are forced to this color
  # (implemented by mapping zeros to a sentinel value beyond zmax)
  field_zero_col = NULL,          # e.g. "darkblue"
  
  # Optional: skip ellipse plotting where an "ocean mask" says ocean
  # Default behavior if TRUE:
  #   - uses ocean_mask_field if provided
  #   - else uses plot_field when draw_field=TRUE
  # Ocean condition is (mask == 0) or NA.
  skip_ocean_ellipses = FALSE,
  ocean_mask_field = NULL,        # matrix on same grid; ocean is defined by == 0
  
  # scaling from params -> ellipse geometry
  kappa2_to_scale = function(kappa2) 1 / sqrt(kappa2),
  base_radius = 0.01,
  max_radius = 0.40,
  
  # ellipse styling
  ellipse_col = "red",
  ellipse_lwd = 0.8,
  ellipse_fill = NA,
  
  # map styling
  map_col = "black",
  map_lwd = 1.2,
  
  # labels and cosmetics
  add_axes = TRUE,
  xlab = "Longitude",
  ylab = "Latitude",
  main = NULL,

   # legend controls
  legend.shrink = 0.8,
  legend.mar = 6,
  legend.lab = "", 
  xaxt = "s", 
  yaxt = "s"    
) {
  stopifnot(all(c("x","y") %in% names(gridList)))
  stopifnot(is.matrix(kappa2), is.matrix(rho), is.matrix(theta))
  stopifnot(all(dim(kappa2) == dim(rho)), all(dim(kappa2) == dim(theta)))
  
  nx <- length(gridList$x)
  ny <- length(gridList$y)
  
  if (!all(dim(kappa2) == c(nx, ny))) {
    stop("Dimensions of kappa2/rho/theta must match c(length(gridList$x), length(gridList$y)).")
  }
  
  # If skipping ocean ellipses, set up the mask (matrix on same grid)
  if (skip_ocean_ellipses) {
    if (is.null(ocean_mask_field)) {
      if (draw_field && !is.null(plot_field)) {
        ocean_mask_field <- plot_field
      } else {
        stop("skip_ocean_ellipses=TRUE requires ocean_mask_field, or draw_field=TRUE with plot_field provided.")
      }
    }
    if (!is.matrix(ocean_mask_field) || !all(dim(ocean_mask_field) == c(nx, ny))) {
      stop("ocean_mask_field must be a matrix with dimensions c(length(gridList$x), length(gridList$y)).")
    }
  }
  
  if (draw_field) {
    if (is.null(plot_field)) stop("draw_field=TRUE requires plot_field (a matrix).")
    if (!is.matrix(plot_field)) stop("plot_field must be a matrix.")
    if (!all(dim(plot_field) == c(nx, ny))) {
      stop("Dimensions of plot_field must match c(length(gridList$x), length(gridList$y)).")
    }
    
    # --- choose zlim ---
    if (is.null(field_zlim)) {
      if (!is.null(field_zero_col)) {
        # compute range excluding exact zeros (so ocean doesn't compress land range)
        z_nonzero <- plot_field[!is.na(plot_field) & plot_field != 0]
        if (length(z_nonzero) == 0) stop("plot_field has no non-zero values to scale colors.")
        field_zlim <- range(z_nonzero)
      } else {
        field_zlim <- range(plot_field, na.rm = TRUE)
      }
    }
    
    # --- build plot inputs ---
    z_plot       <- plot_field
    col_plot     <- field_col
    breaks_plot  <- NULL
    zlim_plot    <- field_zlim
    
    if (!is.null(field_zero_col)) {
      # Map exact zeros to a sentinel value ABOVE the max so they get the last color
      zmin <- field_zlim[1]
      zmax <- field_zlim[2]
      
      ncol <- length(field_col)
      brks <- seq(zmin, zmax, length.out = ncol + 1)
      
      delta <- (zmax - zmin) / ncol
      if (!is.finite(delta) || delta == 0) delta <- 1  # safety
      sentinel <- zmax + 0.5 * delta
      
      # replace exact zeros with sentinel
      z_plot[z_plot == 0] <- sentinel
      
      # extend breaks and colors by one
      breaks_plot <- c(brks, zmax + delta)
      col_plot    <- c(field_col, field_zero_col)
      
      # zlim must include the extended breaks so the sentinel is binned properly
      zlim_plot <- range(breaks_plot)
    }
    
    # IMPORTANT: explicitly call fields::image.plot
    if (!requireNamespace("fields", quietly = TRUE)) {
      stop("Package 'fields' is required for draw_field=TRUE (uses fields::image.plot).")
    }
    
    fields::image.plot(
      x = gridList$x, y = gridList$y, z = z_plot,
      col = col_plot,
      breaks = breaks_plot,
      zlim = zlim_plot,
      xlim = xlim, ylim = ylim,
      xlab = xlab,
      ylab = ylab,
      main = if (!is.null(main)) main else paste0(field, " + ellipses"),
      axes = add_axes,
      legend.shrink = legend.shrink,
      legend.mar = legend.mar,
      legend.args = list(
        text = legend.lab,
        side = 3, las = 1, line = 0.1, 
        cex = 0.8
      ), 
      xaxt = xaxt, 
      yaxt = yaxt
    )
    
  } else {
    plot(
      NA, NA,
      xlim = xlim, ylim = ylim,
      xlab = xlab,
      ylab = ylab,
      main = if (!is.null(main)) main else "Ellipses",
      axes = add_axes
    )
  }
  
  # World map overlay
  map("world", add = TRUE, col = map_col, lwd = map_lwd)
  
  # Subsample every pixel_spacing-th pixel
  ix <- seq(1, nx, by = pixel_spacing)
  iy <- seq(1, ny, by = pixel_spacing)
  
  for (i in ix) {
    for (j in iy) {
      # Optional: skip ocean ellipses based on mask (DEM == 0)
      if (skip_ocean_ellipses) {
        om <- ocean_mask_field[i, j]
        if (is.na(om) || om == 0) next
      }
      
      x0 <- gridList$x[i]
      y0 <- gridList$y[j]
      
      kap <- kappa2[i, j]
      ang <- theta[i, j]
      r   <- rho[i, j]
      
      if (is.na(kap) || is.na(ang) || is.na(r)) next
      
      kap <- max(kap, 1e-12)
      
      s <- kappa2_to_scale(kap)
      s <- min(s, max_radius / base_radius)
      
      b <- base_radius * s
      a <- b * r
      
      a <- min(a, max_radius)
      b <- min(b, max_radius)
      
      poly <- ellipse_xy_helper(x0, y0, a = a, b = b, ang = ang, n = 60)
      
      polygon(
        poly$x, poly$y,
        border = ellipse_col,
        lwd = ellipse_lwd,
        col = ellipse_fill
      )
    }
  }
  
  invisible(list(
    pixel_spacing = pixel_spacing,
    x_idx = ix,
    y_idx = iy
  ))
}


# kappa2 <- params[,,1,chosen_day]
# theta   <- params[,,2,chosen_day]
# rho     <- params[,,3,chosen_day]

# # Classic usage: All ellipses on a white background with map added 

# plot_param_ellipses(
#   kappa2 = kappa2,
#   rho = rho, 
#   theta = theta,
#   gridList = gridList_small, 
#   pixel_spacing = 6,
  
#   draw_field = FALSE, 
#   plot_field = NULL, 
#   field = NULL,
#   field_col = NULL, #topo.colors
  
#   kappa2_to_scale = function(kappa2) 1 / kappa2, #1/ sqrt(kappa2)^1.2
#   base_radius = 0.003, # 0.014
#   max_radius = 0.36,
#   ellipse_col = "red",
#   ellipse_lwd = 1.2,
  
#   map_col = "black",
#   map_lwd = 1.2, 
  
#   main = "Anisotropy ellipses", 
#   xlab = "", 
#   ylab = "", 
#   add_axes = TRUE
#   # xlim = c(10, 14), 
#   # ylim = c(40, 44)
# )




# # Plotting the ellipses on top of a digital elevation map
# # zeroes are assumed to be ocean and set to blue

# plot_param_ellipses(
#   kappa2 = kappa2,
#   rho = rho, 
#   theta = theta,
#   gridList = gridList_small, 
#   pixel_spacing = 6,
  
#   draw_field = TRUE, 
#   plot_field = matrix(cams_df_day$DEM, nrow = 128, ncol = 128), 
#   field = "Elevation",
#   field_col = terrain.colors(256), #topo.colors
#   field_zero_col = "darkblue",
  
#   kappa2_to_scale = function(kappa2) 1 / kappa2, #1/ sqrt(kappa2)^1.2
#   base_radius = 0.003, # 0.014
#   max_radius = 0.36,
#   ellipse_col = "magenta",
#   ellipse_lwd = 1.2,
  
#   map_col = "black",
#   map_lwd = 1.2, 
  
#   main = "Anisotropy ellipses", 
#   xlab = "", 
#   ylab = "", 
#   add_axes = TRUE
# )


# # zooming in on land and removing the huge ocean ellipses
# # this way we can focus on mountains/valleys
# plot_param_ellipses(
#   kappa2 = kappa2,
#   rho = rho, 
#   theta = theta,
#   gridList = gridList_small, 
#   pixel_spacing = 5,
  
#   draw_field = TRUE, 
#   plot_field = matrix(cams_df_day$DEM, nrow = 128, ncol = 128), 
#   field = "Elevation",
#   field_col = terrain.colors(256), #topo.colors
#   field_zero_col = "darkblue",
  
#   skip_ocean_ellipses = TRUE, 
  
#   kappa2_to_scale = function(kappa2) 1 / kappa2, #1/ sqrt(kappa2)^1.2
#   base_radius = 0.02, # 0.014
#   max_radius = 0.4,
#   ellipse_col = "magenta",
#   ellipse_lwd = 1.5,
  
#   map_col = "black",
#   map_lwd = 1.2, 
  
#   main = "Anisotropy ellipses", 
#   xlab = "", 
#   ylab = "", 
#   add_axes = TRUE,
#   xlim = c(6.1, 16),
#   ylim = c(41, 47.9)
# )





# chosen_day <- 103
# average NO2 field and param fields
avg_no2_field <- array(cams_df$NO2, dim = c(128, 128, 365))
avg_no2_field <- apply(avg_no2_field, c(1,2), mean, na.rm = TRUE)

avg_kappa2_field <- apply(params[,,1,], c(1,2), mean, na.rm = TRUE)
avg_theta_field <- apply(params[,,2,], c(1,2), mean, na.rm = TRUE)
avg_rho_field <- apply(params[,,3,], c(1,2), mean, na.rm = TRUE)




# pdf("Paper/ellipses.pdf", width = 9.5, height = 6)
par(mfrow = c(2,3), mar = c(5,4.1,0,6), oma = c(3,0,0,0))

image.plot(
  as.surface(gridList_small, c(avg_no2_field)),
  main = "",
  col = plasma(256), 
  ylab = "Latitude", 
  legend.shrink = 0.8,
  legend.mar = 6,
  legend.args = list(
    tex = expression(mu*g/m^{3}),
    side = 3, 
    las = 1, 
    line = 0.1, 
    cex= 0.8
  ), 
  cex.lab = 1.4, 
  xaxt = "n"
)
world(add = TRUE, col = "white", lwd = 1)


plot_param_ellipses(
  kappa2 = avg_kappa2_field,
  rho = avg_rho_field, 
  theta = avg_theta_field,
  gridList = gridList_small, 
  pixel_spacing = 6,
  
  draw_field = TRUE, 
  plot_field = matrix(cams_df_day$DEM, nrow = 128, ncol = 128), 
  field = "Elevation",
  field_col = terrain.colors(256), #topo.colors
  field_zero_col = "darkblue",
  
  kappa2_to_scale = function(kappa2) 1 / kappa2, #1/ sqrt(kappa2)^1.2
  base_radius = 0.003, # 0.014
  max_radius = 0.36,
  ellipse_col = "magenta",
  ellipse_lwd = 1.2,
  
  map_col = "black",
  map_lwd = 1.2, 
  
  main = "", 
  xlab = "", 
  ylab = "", 
  add_axes = TRUE, 
  legend.lab = "m", 
  xaxt = "n", 
  yaxt = "n"
)


plot_param_ellipses(
  kappa2 = avg_kappa2_field,
  rho = avg_rho_field, 
  theta = avg_theta_field,
  gridList = gridList_small, 
  pixel_spacing = 5,
  
  draw_field = TRUE, 
  plot_field = matrix(cams_df_day$DEM, nrow = 128, ncol = 128), 
  field = "Elevation",
  field_col = terrain.colors(256), #topo.colors
  field_zero_col = "darkblue",
  
  skip_ocean_ellipses = TRUE, 
  
  kappa2_to_scale = function(kappa2) 1 / kappa2, #1/ sqrt(kappa2)^1.2
  base_radius = 0.02, # 0.014
  max_radius = 0.4,
  ellipse_col = "magenta",
  ellipse_lwd = 1.5,
  
  map_col = "black",
  map_lwd = 1.2, 
  
  main = "", 
  xlab = "", 
  ylab = "", 
  add_axes = TRUE,
  xlim = c(6.1, 16),
  ylim = c(41, 47.9), 

  legend.lab = "m", 
  xaxt = "n",
  yaxt = "n"
)


image.plot(
  as.surface(gridList_small, avg_kappa2_field),
  # main = expression(kappa^2),
  col = viridis(256), 
  ylab = "Latitude", xlab = "Longitude",
  legend.shrink = 0.8,
  legend.mar = 6,
  legend.args = list(text = "", side = 3, las = 1, line = 0.1, 
    cex= 0.8), 
  cex.lab = 1.4
)
world(add = TRUE, col = "white", lwd = 1)

image.plot(
  as.surface(gridList_small, avg_theta_field),
  # main = expression(theta), 
  col = viridis(256), 
  xlab = "Longitude",
  legend.shrink = 0.8,
  legend.mar = 6,
  legend.args = list(text = "rad", side = 3, las = 1, line = 0.1, 
    cex= 0.8), 
  cex.lab = 1.4, 
  yaxt = "n"
)
world(add = TRUE, col = "white", lwd = 1)

image.plot(
  as.surface(gridList_small, avg_rho_field),
  # main = expression(rho), 
  col = viridis(256), 
  xlab = "Longitude",
  legend.shrink = 0.8,
  legend.mar = 6,
  legend.args = list(text = "", side = 3, las = 1, line = 0.1, 
    cex= 0.8), 
  cex.lab = 1.4, 
  yaxt = "n"
)
world(add = TRUE, col = "white", lwd = 1)

par(mfrow = c(1,1))
# dev.off()


# ------------------------------------
# CTM Exper Train/Test Figure
# ------------------------------------
chosen_day <- 103

cams_avg <- aggregate(
  NO2 ~ Longitude + Latitude,
  data  = cams_df,
  FUN   = mean,
  na.rm = TRUE
)

stations_avg <- aggregate(
  EEA_NO2 ~ Longitude + Latitude,
  data  = stations_df,
  FUN   = mean,
  na.action = na.pass
)

# subset down to our chosen day 
stations_df_day <- stations_df[
  stations_df$time == (stations_df$time[1] + chosen_day - 1),
]

cams_df_day <- cams_df[
  cams_df$time == (cams_df$time[1] + chosen_day - 1),
]

# scale ssr
cams_df_day[, "ssr"]     <- cams_df_day[, "ssr"] / 10^6
stations_df_day[, "ssr"] <- stations_df_day[, "ssr"] / 10^6

# build sf objects
stations_sf_day <- st_as_sf(
  stations_df_day,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)
cams_sf_day <- st_as_sf(
  cams_df_day,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

stations_sf_avg <- st_as_sf(
  stations_avg,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

cams_sf_avg <- st_as_sf(
  cams_avg,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

nearest_idx <- st_nearest_feature(stations_sf_day, cams_sf_day)
train_idx   <- unique(nearest_idx)

train_df <- cams_df_day[train_idx, ]
test_df  <- cams_df_day[-train_idx, ]

train_df_avg <- cams_avg[train_idx, ]

covars <- c("rh", "ssr", "t2m",
            "windspeed", "sl_blh", "EM_NO2", "lag_cams_no2")

# TRAIN
s <- train_df[, c("Longitude", "Latitude")]
y <- train_df$NO2
Z <- as.matrix(train_df[, covars])

# TEST
s_test <- test_df[, c("Longitude", "Latitude")]
y_test <- test_df$NO2
Z_test <- as.matrix(test_df[, covars])




# actual prediction bit 
# ----- Nonstationary LKinfo: original STUN ns -----
LKinfo_ns <- make_nonstat_LKinfo(
  file_path   = "Italy_AQ_AmortisedLatticeKrig/data/STUN_param_df_2023padded.h5",
  dataset_name= "arx1_surround_30rep_output",
  day         = chosen_day,
  gridlist    = gridList_small,
  normalize   = TRUE,
  sanity_plotting = FALSE,
  sanity_sim      = FALSE
)

# original STUN nonstationary
Lmodel_ns <- LatticeKrig(
  x = s,
  y = y,
  Z = Z,
  LKinfo = LKinfo_ns
)

pred_ns_test <- predict(Lmodel_ns, s_test, Z = Z_test)
pred_ns_full <- predict(
  Lmodel_ns, 
  cams_df_day[, c("Longitude", "Latitude")], 
  as.matrix(cams_df_day[, covars])
)




# plotting setup
nx <- length(gridList_small$x)
ny <- length(gridList_small$y)
n  <- nx * ny

# ---- create full spatial fields with NA elsewhere ----
z_train <- rep(NA_real_, n)
z_test  <- rep(NA_real_, n)
z_test_pred  <- rep(NA_real_, n)

z_train[train_idx] <- y        # actual NO2 values
z_test[-train_idx] <- y_test
z_test_pred[-train_idx] <- pred_ns_test

# ---- reshape to 2D grids ----
zmat_train <- matrix(z_train, nrow = ny, ncol = nx, byrow = FALSE)
zmat_test  <- matrix(z_test,  nrow = ny, ncol = nx, byrow = FALSE)
zmat_test_pred<- matrix(z_test_pred,  nrow = ny, ncol = nx, byrow = FALSE)



zlim <- range(c(y, y_test), na.rm = TRUE)
png("Paper/train_test.png", width = 8, height = 5, units = "in", res = 600)
par(mfrow = c(1,2), mar = c(4,4.1,2,1), bg = "white")
# train points - no legend
image(
  x = gridList_small$x,
  y = gridList_small$y,
  z = zmat_train,
  zlim = zlim,
  col = plasma(256),
  main = "",
  ylab = "Latitude",
  xlab = "Longitude",
  cex.lab = 1.2
)
map("world", add = TRUE, lwd = 1.2)

image(
  x = gridList_small$x,
  y = gridList_small$y,
  z = zmat_test,
  zlim = zlim,
  col = plasma(256),
  main = "",
  xlab = "Longitude",
  ylab = "",
  yaxt = "n",
  cex.lab = 1.2
)
map("world", add = TRUE, lwd = 1.2, col = "white")
par(mfrow = c(1,1), oma = c(0,0,0,0))
dev.off()

png("Paper/train_test_legend.png", width = 5, height = 5, units = "in", res = 600)
plot.new()
# test points - with legend
image.plot(
  x = gridList_small$x,
  y = gridList_small$y,
  z = zmat_test,
  zlim = zlim,
  col = plasma(256),
  main = "",
  xlab = "Longitude",
  ylab = "",
  yaxt = "n",
  legend.shrink = 0.8,
  legend.mar = 6,
  legend.args = list(
    text = expression(mu*g/m^{3}),
    side = 3, las = 1, line = 0.1, cex = 1.0
  ),
  cex.lab = 1.2, 
  legend.only = TRUE
)
dev.off()





# ---- common color scale ----
zlim <- range(c(y, y_test, pred_ns_full), na.rm = TRUE)

par(mfrow = c(1,4), mar = c(4,4.1,2,4), bg = "white")

# FULL field
image.plot(
  as.surface(sGrid_small, cams_df_day$NO2), 
  zlim = zlim,
  col = plasma(256),
  main = "",
  xlab = "Longitude",
  ylab = "Latitude", 
  horizontal = TRUE,
  legend.shrink = 0.8, 
  legend.width = 0.8
)
map("world", add = TRUE, lwd = 1.2, col = "white")

# ---- TRAIN (no colorbar) ----
image.plot(
  x = gridList_small$x,
  y = gridList_small$y,
  z = zmat_train,
  zlim = zlim,
  col = plasma(256),
  main = "Masked Data",
  xlab = "Longitude",
  ylab = "Latitude"
)
map("world", add = TRUE, lwd = 1.2)


# # test 
# image.plot(
#   x = gridList_small$x,
#   y = gridList_small$y,
#   z = zmat_test,
#   zlim = zlim,
#   col = plasma(256),
#   main = "Test",
#   xlab = "Longitude",
#   ylab = "Latitude"
# )
# map("world", add = TRUE, lwd = 1.2, col = "white")

# # predicted test set
# image.plot(
#   x = gridList_small$x,
#   y = gridList_small$y,
#   z = zmat_test_pred,
#   zlim = zlim,
#   col = plasma(256),
#   main = "Predicted Test",
#   xlab = "Longitude",
#   ylab = "Latitude"
# )
# map("world", add = TRUE, lwd = 1.2, col = "white")

# # test set pred differences
# zmat_test_diff <- zmat_test - zmat_test_pred
# max_abs_test_diff <- max(abs(zmat_test_diff), na.rm = TRUE)
# zlim_test_diff <- c(-max_abs_test_diff, max_abs_test_diff)
# 
# image.plot(
#   x = gridList_small$x,
#   y = gridList_small$y,
#   z = zmat_test_diff,
#   zlim = zlim_test_diff,
#   col = cmocean("balance")(256),
#   main = "Predicted Test Diff",
#   xlab = "Longitude",
#   ylab = "Latitude"
# )
# map("world", add = TRUE, lwd = 1.2, col = "black")

# full map pred
image.plot(
  as.surface(gridList_small, pred_ns_full), 
  zlim = zlim,
  col = plasma(256), 
  main = "Predicted NO2 Field", 
  xlab = "Longitude",
  ylab = "Latitude"
)
map("world", add = TRUE, lwd = 1.2, col = "white")


# full map diff
full_diff <- cams_df_day$NO2 - pred_ns_full
max_abs_diff <- max(abs(full_diff), na.rm = TRUE)
zlim_full_diff <- c(-max_abs_diff, max_abs_diff)

image.plot(
  as.surface(gridList_small, full_diff), 
  zlim = zlim_full_diff,
  # col = cmocean("tarn")(256), 
  col = scico(256, palette = 'vik'),
  main = "Difference Field (True - Pred)", 
  xlab = "Longitude",
  ylab = "Latitude"
)
map("world", add = TRUE, lwd = 1.2, col = "black")

par(mfrow = c(1,1))



library(tictoc)
tic()
SE <- predictSE.LKrig(
  Lmodel_ns, 
  cams_df_day[, c("Longitude", "Latitude")], 
  Z = as.matrix(cams_df_day[, covars])
)
toc()

image.plot(as.surface(gridList_small, SE), col = turbo(256))

y_test_LK_def <- cbind(
  y_test_LK_def - 1.96 * sqrt(se_LK_def^2 + LK_def$tau.MLE^2),
  y_test_LK_def,
  y_test_LK_def + 1.96 * sqrt(se_LK_def^2 + LK_def$tau.MLE^2)
)




# Another try 
# -----------------------------
# PLOT CROPPING + AXIS HELPERS
# -----------------------------
xlim_crop <- c(6.4, 18.7)     # <-- adjust as desired
ylim_crop <- c(36.5, 47.1)    # <-- adjust as desired

# for fields::image.plot: suppress colorbar + keep layout stable
no_legend <- list(draw = FALSE)

# y-axis: keep only on first panel
yax_first <- list(side = 2)   # default y-axis
yax_none  <- list(side = 2, labels = FALSE, cex.axis = 0, tck = 0)  # no numbers/ticks

par(
  mfrow = c(1,4), 
  mar = c(3,2.6,2,1.2), 
  oma = c(0,0,2.0,0), 
  mgp   = c(2, 0.6, 0),
  bg = "white"
)

image.plot(
  as.surface(sGrid_small, cams_df_day$NO2),
  zlim = zlim,
  col = plasma(256),
  main = "CTM NO2",
  xlab = "Longitude",
  ylab = "Latitude",
  xlim = xlim_crop,
  ylim = ylim_crop,
  horizontal = TRUE,
  legend.shrink = 0.8, 
  legend.width = 0.8
)

map("world", add = TRUE, lwd = 1.2, col = "white")

image.plot(
  x = gridList_small$x,
  y = gridList_small$y,
  z = zmat_train,
  zlim = zlim,
  col = plasma(256),
  main = "Masked Data",
  xlab = "Longitude",
  ylab = "",
  xlim = xlim_crop,
  ylim = ylim_crop,
  horizontal = TRUE,
  legend.shrink = 0.8, 
  legend.width = 0.8,
  yaxt = "n"
  # axis.args = yax_none
)
# axis(side = 2, labels = TRUE)
map("world", add = TRUE, lwd = 1.2)

image.plot(
  as.surface(gridList_small, pred_ns_full),
  zlim = zlim,
  col = plasma(256),
  main = "Predicted NO2",
  xlab = "Longitude",
  ylab = "",
  xlim = xlim_crop,
  ylim = ylim_crop,
  horizontal = TRUE,
  legend.shrink = 0.8, 
  legend.width = 0.8, 
  yaxt = "n"
)
map("world", add = TRUE, lwd = 1.2, col = "white")

# zlim for diff should reflect only the cropped domain
ix <- gridList_small$x >= xlim_crop[1] & gridList_small$x <= xlim_crop[2]
iy <- gridList_small$y >= ylim_crop[1] & gridList_small$y <= ylim_crop[2]

full_diff_mat <- matrix(full_diff, nrow = ny, ncol = nx, byrow = FALSE)
diff_crop <- full_diff_mat[iy, ix]

max_abs_diff_crop <- max(abs(diff_crop), na.rm = TRUE)
zlim_full_diff_crop <- c(-max_abs_diff_crop, max_abs_diff_crop)

image.plot(
  as.surface(gridList_small, full_diff),
  zlim = zlim_full_diff_crop,
  col  = scico(256, palette = "vik"),
  main = "Difference",
  xlab = "Longitude",
  ylab = "",
  xlim = xlim_crop,
  ylim = ylim_crop,
  horizontal = TRUE,
  legend.shrink = 0.8, 
  legend.width = 0.8, 
  yaxt = "n"
)
map("world", add = TRUE, lwd = 1.2, col = "black")





# One more try  
# -----------------------------

# for fields::image.plot: suppress colorbar + keep layout stable
no_legend <- list(draw = FALSE)

# y-axis: keep only on first panel
yax_first <- list(side = 2)   # default y-axis
yax_none  <- list(side = 2, labels = FALSE, cex.axis = 0, tck = 0)  # no numbers/ticks

pdf("Paper/ctm_reconstruct_exper.pdf", width = 12, height = 5)
par(
  mfrow = c(1,4), 
  mar = c(3,3.6,2,2.2), 
  oma = c(0,0,2.0,0), 
  mgp   = c(2, 0.6, 0),
  bg = "white"
)

image.plot(
  as.surface(sGrid_small, cams_df_day$NO2),
  zlim = zlim,
  col = plasma(256),
  main = "",
  xlab = "Longitude",
  ylab = "Latitude",
  cex.lab = 1.5,
  horizontal = TRUE,
  legend.shrink = 0.7, 
  legend.width = 0.8, 
  legend.args = list(
    text = expression(mu*g/m^{3}),
    side = 4,  las = 1, line = 1, cex = 1.0
  )
)

map("world", add = TRUE, lwd = 1.2, col = "white")

image.plot(
  x = gridList_small$x,
  y = gridList_small$y,
  z = zmat_train,
  zlim = zlim,
  col = plasma(256),
  main = "",
  xlab = "Longitude",
  ylab = "",
  cex.lab = 1.5,
  horizontal = TRUE,
  legend.shrink = 0.7, 
  legend.width = 0.8,
  yaxt = "n",
  legend.args = list(
    text = expression(mu*g/m^{3}),
    side = 4,  las = 1, line = 1, cex = 1.0
  )
  # axis.args = yax_none
)
# axis(side = 2, labels = TRUE)
map("world", add = TRUE, lwd = 1.2)



image.plot(
  as.surface(gridList_small, pred_ns_full),
  zlim = zlim,
  col = plasma(256),
  main = "",
  xlab = "Longitude",
  ylab = "",
  cex.lab = 1.5,
  horizontal = TRUE,
  legend.shrink = 0.7, 
  legend.width = 0.8, 
  yaxt = "n", 
  legend.args = list(
    text = expression(mu*g/m^{3}),
    side = 4,  las = 1, line = 1, cex = 1.0
  )
)
map("world", add = TRUE, lwd = 1.2, col = "white")

# rmse_apr13 <- sqrt(mean((cams_df_day$NO2 - pred_ns_full)^2, na.rm = TRUE))
diff <- cams_df_day$NO2 - pred_ns_full
diffmax <- max(abs(diff), na.rm = TRUE)
zlim_diff <- c(-diffmax, diffmax)

image.plot(
  as.surface(gridList_small, diff),
  zlim = zlim_diff,
  col  = scico(256, palette = "vik"),
  main = "",
  xlab = "Longitude",
  ylab = "",
  cex.lab = 1.5,
  horizontal = TRUE,
  legend.shrink = 0.7, 
  legend.width = 0.8, 
  yaxt = "n", 
  legend.args = list(
    text = expression(mu*g/m^{3}),
    side = 4,  las = 1, line = 1, cex = 1.0
  )
)
map("world", add = TRUE, lwd = 1.2, col = "black")
dev.off()


# ------------------------------------
# CTM Exper Results Figure
# ------------------------------------

results <- readRDS("Italy_AQ_AmortisedLatticeKrig/results/rmse_CTM_reconstruct.rds")
predictions_df <- readRDS("Italy_AQ_AmortisedLatticeKrig/results/predictions_CTM_reconstruct.rds")

# how many days does ns beat stat in terms of RMSE?
sum(results$rmse_test_ns < results$rmse_test_stat, na.rm = TRUE)
sum(results$rmse_test_stat < results$rmse_test_ns, na.rm = TRUE)

total_ns_rmse <- sqrt(mean((predictions_df$ns_pred   - predictions_df$NO2)^2))
total_stat_rmse <- sqrt(mean((predictions_df$stat_pred - predictions_df$NO2)^2))

mean(results$rmse_test_ns, na.rm = TRUE)
mean(results$rmse_test_stat, na.rm = TRUE)

percent_improvement_ns_over_stat <- mean((results$rmse_test_stat - results$rmse_test_ns) / results$rmse_test_stat * 100, na.rm = TRUE)
percent_improvement_ns_over_stat



ns_rmse_field <- aggregate(
  (NO2 - ns_pred)^2 ~ Longitude + Latitude, 
  data = predictions_df, 
  FUN = function(x) sqrt(mean(x, na.rm = TRUE))
)

stat_rmse_field <- aggregate(
  (NO2 - stat_pred)^2 ~ Longitude + Latitude, 
  data = predictions_df, 
  FUN = function(x) sqrt(mean(x, na.rm = TRUE))
)

names(ns_rmse_field)[3] <- "rmse"
names(stat_rmse_field)[3] <- "rmse"

rmse_field_zlim <- range(c(ns_rmse_field$rmse, stat_rmse_field$rmse), na.rm = TRUE)
par(mfrow = c(2,2))
image.plot(
  as.surface(sGrid_small, ns_rmse_field$rmse), 
  col = turbo(256), 
  zlim = rmse_field_zlim
)
image.plot(
  as.surface(sGrid_small, stat_rmse_field$rmse), 
  col = turbo(256), 
  zlim = rmse_field_zlim
)

# set zlim so that the difference plot is centered around zero and has a symmetric color scale
max_abs_diff <- max(abs(ns_rmse_field$rmse - stat_rmse_field$rmse), na.rm = TRUE)
zlim_diff <- c(-max_abs_diff, max_abs_diff)

image.plot(
  as.surface(sGrid_small, ns_rmse_field$rmse - stat_rmse_field$rmse), 
  col = cmocean("balance")(256), 
  zlim = zlim_diff
)

image.plot(
  as.surface(sGrid_small, (stat_rmse_field$rmse - ns_rmse_field$rmse)/stat_rmse_field$rmse), 
  col = cmocean("haline")(256)
)
par(mfrow = c(1,1))



pdf("Paper/ctm_exper_ts_gam.pdf", width = 10, height = 6)
par(mfrow = c(1,2), mar = c(4,5,2,1))

rmse_diff <- results$rmse_test_ns - results$rmse_test_stat

plot(rmse_diff, 
     type = "l", col = "chartreuse3", lwd = 1, 
     main = "", 
     xlab = "Day", 
     ylab = expression(mu*g/m^{3})
     )
points(
  rmse_diff, 
  pch = 21, 
  col = adjustcolor("black", alpha.f = 0.7),     # border color
  bg  = adjustcolor("chartreuse3", alpha.f = 0.6)
)
abline(h = 0, lwd = 3)  # this is the equality line

# # mean line
# abline(
#   h = mean(rmse_diff), 
#   col = "red", 
#   lwd = 3, 
#   lty = 2
# )

# # moving average
# lines(
#   zoo::rollmean(
#     (rmse_diff), 
#     k = 30, 
#     fill = NA, 
#     align = "center"
#   ), 
#   col = "red", 
#   lwd = 2, 
#   lty = 2
# )

# # loess 
day <- seq_along(rmse_diff)
# lo_fit <- loess(rmse_diff ~ day, span = 0.55, degree = 2)  # tweak span as desired
# lo_pred <- predict(lo_fit, newdata = data.frame(day = day))

# lines(lo_pred, col = "red", lwd = 3, lty = 2)

# gam 
gam_fit <- gam(rmse_diff ~ s(day, k = 20, bs = "cc"), method = "ML")
gam_pred <- predict(gam_fit, newdata = data.frame(day = day))

lines(gam_pred, col = "red", lwd = 3, lty = 2)

pct_difference <- (results$rmse_test_ns - results$rmse_test_stat) / results$rmse_test_stat * 100

# rmse_ratio <- results$rmse_test_ns/results$rmse_test_stat

plot(pct_difference, 
     type = "l", col = "cyan3", lwd = 1, 
     ylab = "%", 
     xlab = "Day"
)
points(
  pct_difference, 
  pch = 21, 
  col = adjustcolor("black", alpha.f = 0.7),     # border color
  bg  = adjustcolor("cyan3", alpha.f = 0.6)
)
abline(h = 1, lwd = 3) # this is the equality line\

# LOTS OF LINE OPTIONS BELOW

# # mean line
# abline(
#   h = mean(rmse_ratio), 
#   col = "red", 
#   lwd = 3, 
#   lty = 2
# )

# # moving average
# lines(
#   zoo::rollmean(
#     (rmse_ratio), 
#     k = 30, 
#     fill = NA, 
#     align = "center"
#   ), 
#   col = "red", 
#   lwd = 2, 
#   lty = 2
# )

# # loess 
day <- seq_along(pct_difference)
# lo_fit <- loess(rmse_ratio ~ day, span = 0.55, degree = 2)  # tweak span as desired
# lo_pred <- predict(lo_fit, newdata = data.frame(day = day))

# lines(lo_pred, col = "red", lwd = 3, lty = 2)

# gam 
gam_fit <- gam(pct_difference ~ s(day, k = 20, bs = "cc"), method = "ML")
gam_pred <- predict(gam_fit, newdata = data.frame(day = day))

lines(gam_pred, col = "red", lwd = 3, lty = 2)

# # smooth spline
# ss_fit  <- smooth.spline(x = day, y = rmse_ratio, spar = 1.0)  # spar in [0,1]; smaller = wigglier
# ss_pred <- predict(ss_fit, x = day)$y

# lines(ss_pred, col = "red", lwd = 3, lty = 2)


# # smooth spline
# ss_fit  <- smooth.spline(x = day, y = rmse_diff, spar = 0.7)  # spar in [0,1]; smaller = wigglier
# ss_pred <- predict(ss_fit, x = day)$y

# lines(ss_pred, col = "red", lwd = 3, lty = 2)
par(mfrow = c(1,1))
dev.off()




pdf("Paper/ctm_exper_spatial_rmse.pdf", width = 13, height = 5)
par(mfrow = c(1,3), mar = c(5.6,4,1,6.5), oma = c(0,0,0,0))

image.plot(
  as.surface(sGrid_small, ns_rmse_field$rmse), 
  # make the color rcolorbrewer reds 256
  col = colorRampPalette(brewer.pal(9, "Oranges"))(256),
  # zlim = zlim_diff, 
  main = "", 
  ylab = "Latitude", 
  xlab = "Longitude", 
  cex.lab = 1.5,
  legend.shrink = 0.8, 
  legend.mar = 6,
  legend.args = list(text = expression(mu*g/m^{3}), side = 3, las = 1, line = 0.1, 
    cex= 1.0)
)
world(add = TRUE, col = "black", lwd = 1)

image.plot(
  as.surface(sGrid_small, ns_rmse_field$rmse/avg_no2_field), 
  # make the color rcolorbrewer reds 256
  col = colorRampPalette(brewer.pal(9, "RdPu"))(256),
  # zlim = zlim_diff, 
  main = "", 
  ylab = "", 
  xlab = "Longitude", 
  cex.lab = 1.5,
  legend.shrink = 0.8, 
  legend.mar = 6,
  legend.args = list(text = "", side = 3, las = 1, line = 0.1, 
    cex= 1.0), 
  yaxt = "n"
)
world(add = TRUE, col = "black", lwd = 1)

zlim_diff_max <- max(abs(ns_rmse_field$rmse - stat_rmse_field$rmse), na.rm = TRUE)
zlim_diff <- c(-zlim_diff_max, zlim_diff_max)

#make me a divergent color palette centered at zero where green is negative and red is positive, with 256 colors
# what is the red that rcolorbrewer uses for the highest value in the "Reds" palette? I think it's "#B2182B"
divergent_palette <- colorRampPalette(c("darkgreen", "white", "#B2182B"))(256)

image.plot(
  as.surface(sGrid_small, ns_rmse_field$rmse - stat_rmse_field$rmse), 
  col = divergent_palette, #cmocean("balance")(256)
  zlim = zlim_diff, 
  main = "",
  ylab = "", 
  xlab = "Longitude",
  yaxt = "n",
  cex.lab = 1.5,
  legend.shrink = 0.8, 
  legend.mar = 6,
  legend.args = list(text = expression(mu*g/m^{3}), side = 3, las = 1, line = 0.1, 
    cex= 1.0)
)
world(add = TRUE, col = "black", lwd = 1)

par(mfrow = c(1,1))
dev.off()


# what is the middle color of colorRampPalette(brewer.pal(9, "RdPu"))(256)
color <- colorRampPalette(brewer.pal(9, "RdPu"))(256)[128] #index 128, #ff4690
color

pdf("Paper/ctm_exper_spatiotemporal_rmse.pdf", width = 12, height = 4.8)
par(mfrow = c(1,3), mar = c(5.5,4.2,1,6.4), oma = c(0,0,0,0))


pct_difference <- (results$rmse_test_ns - results$rmse_test_stat) / results$rmse_test_stat * 100

plot(pct_difference, 
     type = "l", 
     col = "grey50", 
     lwd = 1, 
     ylab = "%", 
     xlab = "Day", 
     cex.lab = 1.5
)
points(
  pct_difference, 
  pch = 21, 
  col = adjustcolor("black", alpha.f = 0.7),     # border color
  bg  = adjustcolor("grey50", alpha.f = 0.6)
)
abline(h = 1, lwd = 3) # this is the equality line\

day <- seq_along(pct_difference)

# gam 
gam_fit <- gam(pct_difference ~ s(day, k = 20, bs = "cc"), method = "ML")
gam_pred <- predict(gam_fit, newdata = data.frame(day = day))

lines(gam_pred, col = "#ff4690", lwd = 3)#, lty = 2)


zlim_diff_max <- max(abs(ns_rmse_field$rmse - stat_rmse_field$rmse), na.rm = TRUE)
zlim_diff <- c(-zlim_diff_max, zlim_diff_max)

#make me a divergent color palette centered at zero where green is negative and red is positive, with 256 colors
# what is the red that rcolorbrewer uses for the highest value in the "Reds" palette? I think it's "#B2182B"
divergent_palette <- colorRampPalette(c("darkgreen", "white", "#B2182B"))(256)

image.plot(
  as.surface(sGrid_small, (ns_rmse_field$rmse - stat_rmse_field$rmse)), 
  col = divergent_palette, #cmocean("balance")(256)
  zlim = zlim_diff, 
  main = "",
  ylab = "", 
  xlab = "Longitude",
  yaxt = "n",
  cex.lab = 1.5,
  legend.shrink = 0.8, 
  legend.mar = 6,
  # horizontal = TRUE,
  legend.args = list(text = expression(mu*g/m^{3}), side = 3, las = 1, line = 0.1, 
    cex= 1.0)
)
world(add = TRUE, col = "black", lwd = 1)

image.plot(
  as.surface(sGrid_small, ns_rmse_field$rmse/avg_no2_field), 
  # make the color rcolorbrewer reds 256
  col = colorRampPalette(brewer.pal(9, "RdPu"))(256),
  # zlim = zlim_diff, 
  main = "", 
  ylab = "", 
  xlab = "Longitude", 
  cex.lab = 1.5,
  legend.shrink = 0.8, 
  legend.mar = 6,
  # horizontal = TRUE,
  legend.args = list(text = "", side = 3, las = 1, line = 0.1, 
    cex= 1.0), 
  yaxt = "n"
)
world(add = TRUE, col = "black", lwd = 1)

par(mfrow = c(1,1))
dev.off()



# some more plot options

plot((results$rmse_test_ns - results$rmse_test_stat)/results$rmse_test_stat, 
     type = "l", col = "magenta", lwd = 1, 
     main = "RMSEs % Diff", 
     xlab = "Day", 
     ylab = "RMSE_Nonstat - RMSE_Stat/RMSE_Stat"
)
points(
  (results$rmse_test_ns - results$rmse_test_stat)/results$rmse_test_stat, 
  pch = 21, 
  col = adjustcolor("black", alpha.f = 0.7),     # border color
  bg  = adjustcolor("magenta", alpha.f = 0.6)
)
abline(h = 0, lwd = 3) # this is the equality line
# mean line
abline(
  h = mean((results$rmse_test_ns - results$rmse_test_stat)/results$rmse_test_stat), 
  col = "red", 
  lwd = 3, 
  lty = 2
)




########################################################
############    Station Data            ################
########################################################
stations_df <- readRDS("Italy_AQ_AmortisedLatticeKrig/data/eea_df_2023.rds")
stations_df$ssr <- stations_df$ssr/10^6

stations_df_day <- stations_df[stations_df$time == (19358 + chosen_day - 1), ]

# think im getting white bubbles for NAs 
# but you would know how to fix this better than me 
par(mfrow = c(1,1))
bubblePlot(
  x = stations_df_day$Longitude,
  y = stations_df_day$Latitude,
  z = stations_df_day$EEA_NO2,
  main = paste0("Station NO2 Day ", chosen_day),
  xlab = "Longitude",
  ylab = "Latitude",
  col = plasma(256, alpha = 0.8)
)
world(add = TRUE, col = "black", lwd = 1.2)




########################################################
############    Station CV TS            ###############
########################################################


Y <- read.csv("Italy_AQ_AmortisedLatticeKrig/results/full_crossval_results.csv")
# Y <- read.csv("Italy_AQ_AmortisedLatticeKrig/results/10fold_results.csv")
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


station_metric <- function(df, fit_col, lwr_col, upr_col) {
  df %>%
    mutate(
      sq_err  = (EEA_NO2 - .data[[fit_col]])^2,
      covered = EEA_NO2 >= .data[[lwr_col]] & EEA_NO2 <= .data[[upr_col]],
      piw     = .data[[upr_col]] - .data[[lwr_col]]
    ) %>%
    group_by(AirQualityStation, Longitude, Latitude) %>%
    summarise(
      mean_no2 = mean(EEA_NO2,   na.rm = TRUE),
      rmse = sqrt(mean(sq_err,   na.rm = TRUE)),
      picp = mean(covered,       na.rm = TRUE),
      mpiw = mean(piw,           na.rm = TRUE),
      .groups = "drop"
    )
}

# ── build one wide df: one row per station, metrics for each model ────────────

models <- list(
  lm           = list(fit = "lm_fit",           lwr = "lm_lwr",            upr = "lm_upr"),
  LK_stat      = list(fit = "LK_stat_fit",       lwr = "LK_stat_lwr",       upr = "LK_stat_upr"),
  LK_stun      = list(fit = "LK_stun_fit",       lwr = "LK_stun_lwr",       upr = "LK_stun_upr"),
  LK_stun_adj  = list(fit = "LK_stun_adj_fit",   lwr = "LK_stun_adj_lwr",   upr = "LK_stun_adj_upr")
)

station_wide <- purrr::imap(models, function(cols, name) {
  station_metric(cv_df, cols$fit, cols$lwr, cols$upr) %>%
    rename_with(~ paste0(., "_", name), c(rmse, picp, mpiw, mean_no2))
}) %>%
  purrr::reduce(full_join, by = c("AirQualityStation", "Longitude", "Latitude"))

# collapse the four identical mean_no2_* columns into one
station_wide <- station_wide %>%
  mutate(mean_no2 = rowMeans(select(., starts_with("mean_no2_")), na.rm = TRUE)) %>%
  select(-starts_with("mean_no2_"))


add_alpha <- function(cols, alpha = 0.8) {
  adjustcolor(cols, alpha.f = alpha)
}
colormap1 <- rev(scico(256, palette = "batlow", alpha = 0.8))
colormap2 <- scico(256, palette = "lipari", alpha = 0.8)
colormap3 <- add_alpha(colorRampPalette(brewer.pal(9, "RdPu"))(256), alpha = 0.8)


pdf("Paper/station_bubbleplot.pdf", width = 11, height = 5.6)
par(
  mfrow = c(1,3), 
  mar = c(3,3.6,2,2.2), 
  oma = c(0,0,2.0,0), 
  mgp   = c(2, 0.6, 0),
  bg = "white"
)

bubblePlot(
  x = station_wide$Longitude,
  y = station_wide$Latitude,
  z = station_wide$rmse_LK_stun_adj/station_wide$mean_no2, # relative rmse
  main = "",
  xlab = "Longitude",
  ylab = "Latitude",
  cex.lab = 1.5,
  col = colormap3, 
  horizontal = TRUE,
  legend.shrink = 0.7, 
  legend.width = 0.8, 
  size = 1.5
)
world(add = TRUE)

# nonstat adjusted
bubblePlot(
  x = station_wide$Longitude,
  y = station_wide$Latitude,
  z = station_wide$picp_LK_stun_adj*100,
  main = "",
  xlab = "Longitude",
  ylab = "",
  yaxt = "n",
  cex.lab = 1.5,
  col = colormap1, 
  horizontal = TRUE,
  legend.shrink = 0.7, 
  legend.width = 0.8, 
  legend.args = list(
    text = "%",
    side = 4, las = 1, line = 0.1, cex = 1.0
  ), 
  size = 1.5
)
world(add = TRUE)

bubblePlot(
  x = station_wide$Longitude,
  y = station_wide$Latitude,
  z = station_wide$mpiw_LK_stun_adj,
  main = "",
  xlab = "Longitude",
  ylab = "",
  yaxt = "n",
  cex.lab = 1.5,
  col = colormap2, 
  horizontal = TRUE,
  legend.shrink = 0.7, 
  legend.width = 0.8, 
  size = 1.5
)
world(add = TRUE)
par(mfrow = c(1,1))
dev.off()




bubblePlot(
  x = station_wide$Longitude,
  y = station_wide$Latitude,
  z = (station_wide$rmse_LK_stun_adj - station_wide$rmse_LK_stat)/station_wide$mean_no2, # relative rmse
  main = "Relative RMSE",
  xlab = "Longitude",
  ylab = "Latitude",
  col = turbo(256)
)
world(add = TRUE)

# set green to be values where stun_adj is better, red where stat is better. 
# make it binary, then plot it
diff <- station_wide$rmse_LK_stun_adj - station_wide$rmse_LK_stat
max_abs_diff <- max(abs(diff), na.rm = TRUE)
binary_diff <- ifelse(diff < 0, 1, 0)
colormap_binary <- c("darkgreen", "darkred")

bubblePlot(
  x = station_wide$Longitude,
  y = station_wide$Latitude,
  z = binary_diff,
  main = "LK_stun_adj vs LK_stat",
  xlab = "Longitude",
  ylab = "Latitude",
  col = colormap_binary
)


