#Libraries
library(rhdf5)
library(LatticeKrig)
library(tictoc)

setwd("C:/Users/anton/Desktop/Research/Italy_AQ_paper/Italy_AQ_AmortisedLatticeKrig")

#---------------------------------
# Fitting ARX1 and LM to data 
#---------------------------------
# cams_df_old <- readRDS("data/cams_df_2023.rds")
# cams_df_new <- readRDS("data/cams_df_2023padded.rds")
# cams_df_new_sub <- cams_df_new[
#   cams_df_new$time >= 19358 & cams_df_new$time <= 19722, 
# ]
# test.for.zero(cams_df_old$NO2, cams_df_new_sub$NO2)
# test.for.zero(cams_df_old$t2m, cams_df_new_sub$t2m)
# test.for.zero(cams_df_old$DEM, cams_df_new_sub$DEM)
# 
# gridList_small <- list( 
#   x = seq( 6.1,18.8,by=.1), 
#   y= seq( 35.2,47.9,by=.1)
# )
# 
# par(mfrow = c(1,3))
# imagePlot(as.surface(gridList_small, cams_df_old$DEM[cams_df_old$time == 19460]), 
#           col=cividis(256), 
#           main = "old"  ) 
# imagePlot(as.surface(gridList_small, cams_df_new_sub$DEM[cams_df_new_sub$time == 19460]), 
#           col=cividis(256), 
#           main = "new"  ) 
# imagePlot(
#   as.surface(
#     gridList_small, 
#     cams_df_new_sub$DEM[cams_df_new_sub$time == 19460] 
#     - cams_df_old$DEM[cams_df_old$time == 19460]
#     ), 
#   col=turbo(256), 
#   main = "diff"  
# ) 
# par(mfrow = c(1,1))

cams_df <- readRDS("data/cams_df_2023padded.rds")
# first day does not have lagged NO2 (19336)
cams_df <- cams_df[cams_df$time >= 19337, ] 

#split by day 
cams_by_day <- split(cams_df, cams_df$time)

# fit an ARX(1) model for each day
arx1_by_day <- lapply(
  cams_by_day,
  function(df_day) {
    lm(
      NO2 ~ Longitude + Latitude + rh + ssr + t2m + windspeed + sl_blh +
        EM_NO2 + DEM + lag_cams_no2, data = df_day)
  }
)

# #fit lm model for each day 
# lm_by_day <- lapply(
#   cams_by_day,
#   function(df_day) {
#     lm(
#       NO2 ~ Longitude + Latitude + rh + ssr + t2m + windspeed + sl_blh +
#         EM_NO2, data = df_day)
#   }
# )

rm(cams_by_day)

# add results to cams df 
# cams_df$lm_predict <- NA
# cams_df$lm_resid <- NA
cams_df$arx1_predict <- NA
cams_df$arx1_resid <- NA

#coefs to compare to eea model
arx1_coefs <- data.frame(
  Intercept = NA, 
  Longitude = NA,  
  Latitude  = NA, 
  rh  = NA, 
  ssr = NA, 
  t2m = NA, 
  windspeed = NA, 
  sl_blh = NA, 
  EM_NO2 = NA, 
  DEM = NA, 
  lag_cams_no2 = NA
)

# lm_coefs <- data.frame(
#   Intercept = NA, 
#   Longitude = NA,  
#   Latitude  = NA, 
#   rh  = NA, 
#   ssr = NA, 
#   t2m = NA, 
#   windspeed = NA, 
#   sl_blh = NA, 
#   EM_NO2 = NA
# )

# i <- 1
# for (tt in names(arx1_by_day)) {
#   # match the df 
#   idx <- cams_df$time == as.numeric(tt)
#   
#   # cams_df$lm_predict[idx] <- lm_by_day[[i]]$fitted.values
#   # cams_df$lm_resid[idx] <- lm_by_day[[i]]$residuals
#   
#   cams_df$arx1_predict[idx] <- arx1_by_day[[i]]$fitted.values
#   cams_df$arx1_resid[idx] <- arx1_by_day[[i]]$residuals
#   
#   # coef dfs
#   arx1_coefs[i, ] <- arx1_by_day[[i]]$coefficients
#   # lm_coefs[i, ] <- lm_by_day[[i]]$coefficients
#   i <- i + 1
# }

i <- 1
for (tt in names(arx1_by_day)) {
  
  # convert the loop key to Date
  tt_date <- as.Date(tt)
  
  # match the df (make sure cams_df$time is Date too)
  idx <- as.Date(cams_df$time) == tt_date
  
  cams_df$arx1_predict[idx] <- arx1_by_day[[i]]$fitted.values
  cams_df$arx1_resid[idx]   <- arx1_by_day[[i]]$residuals
  
  arx1_coefs[i, ] <- arx1_by_day[[i]]$coefficients
  i <- i + 1
}



#---------------------------------
# Checking results with plotting
#---------------------------------

rows <- 128
cols <- 128

gridList <- list( 
  x = seq( 6.1,18.8,by=.1), 
  y= seq( 35.2,47.9,by=.1)
)
sGrid<- make.surface.grid(gridList)


# choose a time
time_idx <- 19358
time_idx <- 19460


# first, predictions
par(mfrow = c(1,2))
pred_max <- max(
  cams_df$NO2[cams_df$time == time_idx], 
  # cams_df$lm_predict[cams_df$time == time_idx], 
  cams_df$arx1_predict[cams_df$time == time_idx]
)
pred_min <- min(
  cams_df$NO2[cams_df$time == time_idx], 
  # cams_df$lm_predict[cams_df$time == time_idx], 
  cams_df$arx1_predict[cams_df$time == time_idx]
)
imagePlot(as.surface(gridList, cams_df$NO2[cams_df$time == time_idx]), col=plasma(256), 
          main = "Original",  zlim = c(pred_min, pred_max)) 

# imagePlot(as.surface(gridList, cams_df$lm_predict[cams_df$time == time_idx]), col=plasma(256), 
#           main = "lm predict", zlim = c(pred_min, pred_max)) 

imagePlot(as.surface(gridList, cams_df$arx1_predict[cams_df$time == time_idx]), col=plasma(256), 
          main = "arx1 predict", zlim = c(pred_min, pred_max)) 
par(mfrow = c(1,1))


# now some residuals
par(mfrow = c(1,1))
res_max <- max(
  # cams_df$lm_resid[cams_df$time == time_idx], 
  cams_df$arx1_resid[cams_df$time == time_idx]
)
res_min <- min(
  # cams_df$lm_resid[cams_df$time == time_idx], 
  cams_df$arx1_resid[cams_df$time == time_idx]
)

# imagePlot(as.surface(gridList, cams_df$lm_resid[cams_df$time == time_idx]), col=turbo(256), 
#           main = "lm resid")

imagePlot(as.surface(gridList, cams_df$arx1_resid[cams_df$time == time_idx]), col=turbo(256), 
          main = "arx1 resid")#, zlim = c(res_min, res_max)) 
par(mfrow = c(1,1))

# final save of data
saveRDS(cams_df, file = "data/cams_df_2023padded.rds")
# 21 extra days in the beginning
# 15 extra days in the end

#---------------------------------
# saving residuals as hdf5 files 
#---------------------------------

#filepath to where the h5 file will be saved 
h5_file_path <- "data/NO2_resid_df_2023padded.h5"

# create the file
h5createFile(h5_file_path)

rows <- sqrt(nrow(cams_df)/length(unique(cams_df$time)))
cols <- rows

# lm_resid_forh5 <- array(
#   data = cams_df$lm_resid, 
#   dim = c(rows, cols, length(unique(cams_df$time)))
# )

arx1_resid_forh5 <- array(
  data = cams_df$arx1_resid, 
  dim = c(rows, cols, length(unique(cams_df$time)))
)

# h5write(obj = lm_resid_forh5, file = h5_file_path, name = "lm_resid")  
h5write(obj = arx1_resid_forh5, file = h5_file_path, "arx1_resid") 

h5ls(h5_file_path)

H5close()



# OBSOLETE AFTER THIS, CHECKING THAT DEM MADE A DIFF

# old_cams_df <- readRDS("data/old/cams_df_2023.rds")
# cams_df_subset <- cams_df[
#   cams_df$time >= 19358 & cams_df$time <= 19722, 
# ]
# 
# test.for.zero(old_cams_df$arx1_resid, cams_df_subset$arx1_resid)
# 
# par(mfrow = c(1,3))
# imagePlot(as.surface(gridList, old_cams_df$arx1_resid[old_cams_df$time == 19434]), 
#           col=turbo(256), 
#           main = "old"  )
# imagePlot(as.surface(gridList, cams_df_subset$arx1_resid[cams_df_subset$time == 19434]), 
#           col=turbo(256), 
#           main = "new"  )
# # plot their diff
# imagePlot(as.surface(gridList, 
#           old_cams_df$arx1_resid[old_cams_df$time == 19434] - 
#             cams_df_subset$arx1_resid[cams_df_subset$time == 19434]), 
#           col=turbo(256), 
#           main = "diff"  )
# par(mfrow = c(1,1))

# 
# # BEYOND THIS IS OBSOLETE
# 
# #---------------------------------
# # GLS Fit 
# #---------------------------------
# 
# # crack open updated cams_df 
# cams_df <- readRDS("data/cams_df_2023.rds")
# 
# cams_by_day <- split(cams_df, cams_df$time)
# 
# # need this for nonstationary awghts 
# predict.multivariateSurfaceGrid<- function(object,x){
#   dimZ<- dim( object$z)
#   L<- dimZ[3]
#   out<- matrix( NA, nrow= nrow(x), ncol=L)
#   for (  l in 1:L){
#     out[,l]<- interp.surface( 
#       list( x=object$x,y=object$y, z=object$z[,,l]) , x)
#   }
#   return( out)
# }
# 
# # function to make nonstat lkinfo 
# make_nonstat_LKinfo <- function(
#     file_path, 
#     dataset_name,
#     day, 
#     gridlist,
#     normalize = TRUE, 
#     NC = 128, 
#     nlevel = 1, 
#     NC.buffer = 0, 
#     sanity_plotting = FALSE, 
#     sanity_sim = FALSE
# ){
#   params <- h5read(file_path, dataset_name)
#   H5close()
#   # for plotting and image purposes, need to flip upside down
#   params <- params[,128:1,,day]
#   
#   # recover kappa
#   kappa2 <- exp(params[,,1])
#   # just in case we need awght
#   awght <- kappa2 + 4
#   # theta needs to be transformed like this 
#   theta <- params[,,2] + pi/2
#   rho <- params[,,3]
#   
#   if (sanity_plotting){
#     par(mfrow = c(1,3))
#     imagePlot(as.surface(gridlist, kappa2), main="kappa2", col = viridis(256))
#     world(add=TRUE, col = "white", lwd = 1)
#     imagePlot(as.surface(gridlist, theta-pi/2), main="theta (adj)", col = viridis(256))
#     world(add=TRUE, col = "white", lwd = 1)
#     imagePlot(as.surface(gridlist, rho), main="rho", col = viridis(256))
#     world(add=TRUE, col = "white", lwd = 1)
#     par(mfrow = c(1,1))
#   }
#   
#   # need these for encoding into LK
#   rhox <- sqrt(rho)
#   rhoy <- 1/rhox
#   
#   # create H tensor out of params
#   H11 <- ( rhox^2 * (cos(theta))^2) + ( rhoy^2 * (sin(theta))^2 ) 
#   H12 <- (rhoy^2 - rhox^2)*(sin(theta)*cos(theta))
#   H21 <- H12 
#   H22 <- (rhox^2 * (sin(theta))^2) + (rhoy^2 * (cos(theta))^2)
#   
#   rows <- length(gridlist$x)
#   
#   # fill the high dimensional stencil (9 fields)
#   stencil_tensor <- array( NA, c( rows,rows,9))
#   stencil_tensor[,,1] <- 0.5 * H12
#   stencil_tensor[,,2] <- -H22
#   stencil_tensor[,,3] <- -0.5 * H12
#   stencil_tensor[,,4] <- -H11
#   stencil_tensor[,,5] <- kappa2 + 2 * H11 + 2 * H22
#   stencil_tensor[,,6] <- -H11
#   stencil_tensor[,,7] <- -0.5 * H12
#   stencil_tensor[,,8] <- -H22
#   stencil_tensor[,,9] <- 0.5 * H12
#   
#   # next, we put everything into awght obj of a particular class
#   awght_obj <- list( x= gridlist$x,  y= gridlist$y, z=stencil_tensor )
#   class( awght_obj)<- "multivariateSurfaceGrid"
#   
#   sGrid <- make.surface.grid(gridlist)
#   
#   LKinfo <- LKrigSetup(
#     # dont change grid
#     sGrid, 
#     # you change awght indirectly with day and file selection
#     a.wghtObject =  awght_obj,
#     # the rest of these you change directly in function calls 
#     NC = NC, 
#     nlevel = nlevel,
#     normalize=normalize,
#     NC.buffer = NC.buffer
#   )
#   
#   if(sanity_sim){
#     test <- LKrig.sim(
#       sGrid,
#       LKinfo=LKinfo,
#       M = 1
#     )
#     
#     if (sanity_plotting){
#       imagePlot(as.surface(gridlist, test), col = turbo(256))
#       world(add=TRUE, col = "black", lwd = 1)
#     }
#   }
#   
#   return(LKinfo)
# }
# 
# rows <- 128
# cols <- 128
# 
# gridList <- list( 
#   x = seq( 6.1,18.8,by=.1), 
#   y= seq( 35.2,47.9,by=.1)
# )
# sGrid<- make.surface.grid(gridList)
# 
# # formula for the gls 
# gls_formula <- NO2 ~ Longitude + Latitude + rh + ssr + t2m + windspeed + 
#   sl_blh + EM_NO2 + lag_cams_no2
# 
# # time vals to loop through 
# time_vals <- sort(unique(cams_df$time))  # should be length 365
# 
# cams_df$arx1_predict4 <- NA_real_
# cams_df$arx1_resid4 <- NA_real_
# 
# # loop through days to do FGLS
# for (day_idx in seq_along(time_vals)) {
#   tic()
#   tt <- time_vals[day_idx]  # actual time value in cams_df
#   
#   # subset CAMS for this day
#   df_day <- subset(cams_df, time == tt)
#   df_day$ssr <- df_day$ssr/10^6
#   
#   # build nonstationary LKinfo for this day
#   LKinfo_ns <- make_nonstat_LKinfo(
#     file_path       = "data/STUN_param3_df_2023.h5",
#     dataset_name    = "arx1_surround3_30rep_output",
#     day             = day_idx,        # NOTE: indexes 1..365
#     gridlist        = gridList,
#     normalize       = TRUE,
#     sanity_plotting = FALSE,
#     sanity_sim      = FALSE
#   )
#   
#   # precision matrix for the spatial field (spam object)
#   Q <- LKrig.precision(LKinfo_ns)  # dimension 16384 x 16384
#   
#   # design matrix X and response y for this day
#   X <- model.matrix(gls_formula, data = df_day)
#   y <- df_day$NO2
#   
#   # GLS estimate: beta_hat = (X^T Q X)^(-1) X^T Q y
#   # Use spam-friendly operations
#   # Q %*% X gives an n x p dense matrix (n ~ 16384, p ~ 10)
#   QX   <- Q %*% X                # n x p
#   Qy   <- Q %*% y                # n x 1
#   XtQX <- crossprod(X, QX)       # p x p
#   XtQy <- crossprod(X, Qy)       # p x 1
#   
#   beta_gls <- solve(XtQX, XtQy)  # p x 1
#   
#   # fitted values and residuals
#   y_hat <- as.numeric(X %*% beta_gls)
#   resid <- y - y_hat
#   
#   # write back into cams_df for the rows of this day
#   idx <- which(cams_df$time == tt)
#   cams_df$arx1_predict4[idx] <- y_hat
#   cams_df$arx1_resid4[idx]   <- resid
#   
#   cat("Finished FGLS for day index", day_idx, "time", tt, "\n")
#   toc()
# }
# 
# 
# 
# par(mfrow = c(3,4))
# par(mar = c(2,2,2,2))
# chosen_day <- 19358
# 
# plotmax <- max(
#   cams_df$arx1_resid[cams_df$time == chosen_day], 
#   cams_df$arx1_resid1[cams_df$time == chosen_day], 
#   cams_df$arx1_resid2[cams_df$time == chosen_day], 
#   cams_df$arx1_resid3[cams_df$time == chosen_day], 
#   cams_df$arx1_resid4[cams_df$time == chosen_day]
# )
# plotmin <- min(
#   cams_df$arx1_resid[cams_df$time == chosen_day], 
#   cams_df$arx1_resid1[cams_df$time == chosen_day], 
#   cams_df$arx1_resid2[cams_df$time == chosen_day],
#   cams_df$arx1_resid3[cams_df$time == chosen_day], 
#   cams_df$arx1_resid4[cams_df$time == chosen_day]
# )
# imagePlot(as.surface(gridList, cams_df$arx1_resid[cams_df$time == chosen_day]), 
#           col=turbo(256), main = "arx1 resid", zlim = c(plotmin, plotmax))
# imagePlot(as.surface(gridList, cams_df$arx1_resid1[cams_df$time == chosen_day]), 
#           col=turbo(256), main = "arx1 gls resid 1", zlim = c(plotmin, plotmax))
# imagePlot(as.surface(gridList, 
#           cams_df$arx1_resid1[cams_df$time == chosen_day] - 
#             cams_df$arx1_resid[cams_df$time == chosen_day]), 
#           col=turbo(256), main = "diff1")
# 
# imagePlot(as.surface(gridList, cams_df$arx1_resid1[cams_df$time == chosen_day]), 
#           col=turbo(256), main = "arx1 gls resid 1", zlim = c(plotmin, plotmax))
# imagePlot(as.surface(gridList, cams_df$arx1_resid2[cams_df$time == chosen_day]), 
#           col=turbo(256), main = "arx1 gls resid 2", zlim = c(plotmin, plotmax))
# imagePlot(as.surface(gridList, 
#           cams_df$arx1_resid2[cams_df$time == chosen_day] - 
#             cams_df$arx1_resid1[cams_df$time == chosen_day]), 
#           col=turbo(256), main = "diff2")
# 
# 
# imagePlot(as.surface(gridList, cams_df$arx1_resid2[cams_df$time == chosen_day]), 
#           col=turbo(256), main = "arx1 gls resid 2", zlim = c(plotmin, plotmax))
# imagePlot(as.surface(gridList, cams_df$arx1_resid3[cams_df$time == chosen_day]), 
#           col=turbo(256), main = "arx1 gls resid 3", zlim = c(plotmin, plotmax))
# imagePlot(as.surface(gridList, 
#           cams_df$arx1_resid3[cams_df$time == chosen_day] - 
#             cams_df$arx1_resid2[cams_df$time == chosen_day]), 
#           col=turbo(256), main = "diff3")
# 
# imagePlot(as.surface(gridList, cams_df$arx1_resid3[cams_df$time == chosen_day]), 
#           col=turbo(256), main = "arx1 gls resid 3", zlim = c(plotmin, plotmax))
# imagePlot(as.surface(gridList, cams_df$arx1_resid4[cams_df$time == chosen_day]), 
#           col=turbo(256), main = "arx1 gls resid 4", zlim = c(plotmin, plotmax))
# imagePlot(as.surface(gridList, 
#           cams_df$arx1_resid4[cams_df$time == chosen_day] - 
#             cams_df$arx1_resid3[cams_df$time == chosen_day]), 
#           col=turbo(256), main = "diff4")
# 
# par(mfrow = c(1,1))
# 
# # save updated cams_df
# saveRDS(cams_df, file = "data/cams_df_2023.rds")
# 
# 
# # saving residuals as h5 
# h5_file_path <- "data/NO2_resid4_df_2023.h5"
# 
# # create the file
# h5createFile(h5_file_path)
# 
# rows <- sqrt(nrow(cams_df)/365)
# cols <- rows
# 
# 
# arx1_resid_forh5 <- array(
#   data = cams_df$arx1_resid4, 
#   dim = c(rows, cols, 365)
# )
# 
# h5write(obj = arx1_resid_forh5, file = h5_file_path, "arx1_resid4") 
# 
# h5ls(h5_file_path)
# 
# H5close()
# 
