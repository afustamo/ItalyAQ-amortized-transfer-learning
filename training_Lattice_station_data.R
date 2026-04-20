# try2 is the baseline
# try3  just background stations and concentrations less than 80, Ndata > 250
# try4 is just CV with log scale
# try5 - TO BE DONE - try other kappa2_weights, in particular, use the shapefile of Italy + 20km of buffer!
# or a kind of a smooth, like 1 at coast and gradually toward 0 going far from it
# and we keep also the log scale for the response variable
# try6 - same weights as try5 but without log scale

HPC <- as.logical(Sys.getenv("AQ_amortised_HPC"))
if (is.na(HPC)) {
  HPC <- F
}

# Packages and source codes ####
library(doParallel)
library(foreach)
if (HPC) {
  registerDoParallel(cores = detectCores() - 1)
} else{
  registerDoParallel(cores = detectCores() / 2)
}
library(rhdf5)
library(LatticeKrig)
library(sf)
# library(FRK)
library(sp)
library(gstat)
library(tictoc)

# %% settings

if (HPC) {
  surfaces <- T 
  only_params <- F
  CV_scheme <- T
  total_time <- 365
  try(setwd("AQ_amortised"))
} else{
  surfaces <- T
  CV_scheme <- F
  total_time <- 60
}

cams_df <- readRDS("data/cams_df_2023.rds")
elev_map <-cams_df$DEM[cams_df$time == 19358]
elev_map_mask <- elev_map 
elev_map_mask[elev_map_mask != 0] <- 1
kappa2_weights <- elev_map_mask


# # try5 - try6 without log
# load("data/Italy_shp.Rdata")

# plot(Italy)
# grid_cams <- unique(cams_df[,c("Longitude","Latitude")])
# coordinates(grid_cams) <- c("Longitude","Latitude")
# Italy <- st_as_sf(Italy)
# grid_cams <- st_as_sf(grid_cams)
# st_crs(Italy)
# st_crs(grid_cams) <- st_crs(4326)
# Italy <- st_transform(Italy,st_crs(4326))
# grid_dist <- st_distance(grid_cams, Italy)
# # image(matrix(grid_dist,128,128))
# dist_sp <- SpatialPointsDataFrame(as_Spatial(grid_cams),data.frame(dist_coast=grid_dist))
# dist_sp@data$kappa2_old <- kappa2_weights
# summary(dist_sp@data$dist_coast)
# dist_sp@data$kappa2_weights <- NA
# dist_sp@data$kappa2_weights[dist_sp@data$kappa2_old==0] <- 1
# dist_sp@data$kappa2_weights <- exp(-as.numeric(dist_sp@data$dist_coast/30000))
# dist_sp@data$kappa2_weights[as.numeric(dist_sp@data$dist_coast)>30000] <- 0
# # gridded(dist_sp)<-T
# # spplot(dist_sp,"kappa2_weights")
# kappa2_weights<- dist_sp@data$kappa2_weights
# # ggplot(cams_df[cams_df$time == 19358,],aes(x=Longitude,y=Latitude))+
# #   geom_tile(fill=as.factor(kappa2_weights))


rm(cams_df,grid_cams,dist_sp)

if (CV_scheme) {
  CV_type <- Sys.getenv("AQ_amortised_CV_type")
  if (CV_type == "") {
    CV_type <- "10-fold"
  }
}

# 1. "10-fold"
# look inside "paste0("data/output/", CV_type" , find the biggest try1 try2 etc. and make a new folder for the new try

if (CV_scheme) {
  try_ex <- list.files(paste0("data/output/", CV_type),pattern = "try")
  n_try_ex <- as.numeric(substr(try_ex,4,nchar(try_ex)))
  n_try_ex <- max(n_try_ex)
  n_try_new <- n_try_ex + 1
  n_try_new_cv <- n_try_new
}

if (surfaces) {
  try_ex1 <- list.files(paste0("data/output/map"),pattern = "try")
  try_ex2 <- list.files(paste0("data/output/coef"),pattern = "try")
  n_try_ex1 <- max(as.numeric(substr(try_ex1,4,nchar(try_ex1))))
  n_try_ex2 <- max(as.numeric(substr(try_ex2,4,nchar(try_ex2))))
  n_try_ex <- max(n_try_ex1, n_try_ex2)
  n_try_new <- n_try_ex + 1
  n_try_new_surf <- n_try_new
}

if (surfaces & CV_scheme){
  n_try_new <- max(n_try_new_cv,n_try_new_surf)
}

if (CV_scheme) {
dir.create(paste0("data/output/", CV_type, "/try", n_try_new))
CV_folder <- file.path("data/output", CV_type, paste0("try", n_try_new))
}
if (surfaces) {
dir.create(paste0("data/output/map/try", n_try_new))
dir.create(paste0("data/output/coef/try", n_try_new))
map_folder <- file.path(paste0("data/output/map/try", n_try_new))
coef_folder <- file.path(paste0("data/output/coef/try", n_try_new))
}

if (CV_scheme){
print("CV_folder")
print(CV_folder)
  }
if (surfaces){
print("map_folder")
print(map_folder)
print("coef_folder")
print(coef_folder)
  }

# for HPC
source("R_scripts/hacked_LatticeKrig.R")
source("R_scripts/LK_functions.R")

# Loading datasets ####
stations_df <- readRDS("data/eea_df_2023.rds")
pred_df <- readRDS("data/pred_df_2023.rds")
params_all_surr30 <- h5read("data/STUN_param_df_2023.h5", "arx1_surround_30rep_output")
params_all_surr1 <- h5read("data/STUN_param_df_2023.h5", "arx1_surround_1rep_output")
params_all_rscale30 <- h5read("data/STUN_param_df_2023.h5", "arx1_rscale_30rep_output")
params_all_rscale1 <- h5read("data/STUN_param_df_2023.h5", "arx1_rscale_1rep_output")


# h5ls("data/STUN_param_df_2023.h5")
H5close()

load("data/Station_registry_information.rda")
stations_df <- merge(stations_df,Station_registry_information[,c("AirQualityStation",
"AirQualityStationType")])

#try3 
# change also the number of stations necessary to 250 instead of 450
# hist(stations_df$EEA_NO2,breaks = 150) #try3
# hist(log(stations_df$EEA_NO2),breaks = 150) #try3
# quantile(stations_df$EEA_NO2,.999,na.rm=T)
stations_df$EEA_NO2[stations_df$EEA_NO2 >= 80] <- NA
stations_df <- subset(stations_df,AirQualityStationType=="background")

# #try4 and try5
# stations_df$EEA_NO2 <- log(stations_df$EEA_NO2)

# days involved
station_time <- unique(stations_df$time)
station_time <- station_time[order(station_time)]
aa<-c()
time_sel <- c()
all_N_data <- c()
for (i in 1:365) {
  stations_df_day <- stations_df[stations_df$time == station_time[i], ]
  stations_df_day <- stations_df_day[!is.na(stations_df_day$EEA_NO2), ]
  N_data <- nrow(stations_df_day)
  all_N_data <- c(all_N_data,N_data)
  if(N_data>=250){aa <- c(aa,TRUE)
  time_sel <- c(time_sel,station_time[i])}
}
head(time_sel - 19357)
which(station_time==as.Date("2023-01-10"))

# %%

# Begin cycle ####
foreach(
  time_cycle = 1:total_time,
  .packages = c("rhdf5", "LatticeKrig", "sf", "FRK", "sp", "gstat"),
  .export = c("nonstat_Q","kappa2_weights"),
  .combine = rbind
) %dopar% {
  ## functions ####
  source("R_scripts/hacked_LatticeKrig.R")
  source("R_scripts/LK_functions.R")
  
  predict.multivariateSurfaceGrid <- function(object, x) {
    dimZ <- dim(object$z)
    L <- dimZ[3]
    out <- matrix(NA, nrow = nrow(x), ncol = L)
    for (l in 1:L) {
      out[, l] <- interp.surface(list(
        x = object$x,
        y = object$y,
        z = object$z[, , l]
      ), x)
    }
    return(out)
  }
  # est_additional_kappa2 <- function(s, y, Z, LKinfo, kappa2_weights, verbose = FALSE) {
  #   kappa2_weights <- kappa2_weights
  #   
  #   LKinfo_additional <- LKinfo
  #   LKinfo_additional$a.wght <- NA
  #   LKinfo_additional$a.wghtObject <- NULL
  #   LKinfo_additional$a.wght <- list()
  #   LKinfo_additional$a.wght[[1]] <- 4.01
  #   attr(LKinfo_additional$a.wght, "isotropic") <- TRUE
  #   
  #   Lmodel_additional <- LatticeKrig_hacked(s, y, Z, LKinfo = LKinfo_additional, findAwght = T)
  #   
  #   LKinfo_new <- LKinfo
  #   LKinfo_new$a.wght[[1]][, 5] <- LKinfo_new$a.wght[[1]][, 5] +
  #     (Lmodel_additional$LKinfo$a.wght[[1]] - 4)
  #   
  #   return(LKinfo_new)
  # }
  
  est_additional_kappa2 <- function(
    s, 
    y, 
    Z, 
    LKinfo, 
    kappa2_weights,
    verbose = FALSE
){
  
  kappa2_weights <- kappa2_weights
  # assign("kappa2_weights", as.numeric(kappa2_weights), envir = .GlobalEnv)
  
  LKinfo_additional <- LKinfo
  LKinfo_additional$a.wght <- NA
  LKinfo_additional$a.wghtObject <- NULL
  LKinfo_additional$a.wght <- list()
  LKinfo_additional$a.wght[[1]] <- 4.01
  attr(LKinfo_additional$a.wght, "isotropic") <- TRUE
  
  Lmodel_additional <- LatticeKrig_hacked(
    s,
    y,
    Z,
    LKinfo = LKinfo_additional, 
    findAwght = T
  )
  
  LKinfo_new <- LKinfo
  LKinfo_new$a.wght[[1]][,5] <- LKinfo_new$a.wght[[1]][,5] + 
    (kappa2_weights * (Lmodel_additional$LKinfo$a.wght[[1]]-4))
  
  return(LKinfo_new)
}
  # now let's do the additional kappa2 modification
  LatticeKrig_hacked <- function(x,
                                 y,
                                 Z = NULL,
                                 weights = NULL,
                                 nlevel = 3,
                                 findAwght = FALSE,
                                 LKinfo = NULL,
                                 X = NULL,
                                 U = NULL,
                                 na.rm = TRUE,
                                 tol = .005,
                                 verbose = FALSE,
                                 ...) {
    # a crisp wrapper where many default values are exercised.
    x <- as.matrix(x)
    y <- as.matrix(y)
    if (is.null(weights)) {
      weights <- rep(1, nrow(y))
    }
    # adjust for missing values
    ind <- is.na(y)
    if (any(ind)) {
      if (na.rm) {
        x <- x[!ind, ]
        y <- y[!ind]
        weights <- weights[!ind]
        warning("NAs removed")
        if (!is.null(Z)) {
          Z <- as.matrix(Z)[!ind, ]
        }
      }
      else{
        stop("NAs in y")
      }
    }
    #
    if (is.null(LKinfo)) {
      argList <- list(...)
      # determine the geometry/dimension if not specified
      # set up some thin plate spline like default models for just Euclidean spatial domains
      # in 1,2 and 3 dimensions.
      argList <- LatticeKrigEasyDefaults(argList, nlevel, x)
      if (verbose) {
        cat("extra args:", fill = TRUE)
        print(names(argList))
      }
      LKinfo <- do.call("LKrigSetup", c(list(
        x = x,
        nlevel = nlevel,
        verbose = FALSE
      ), argList))
    }
    if (verbose) {
      print(LKinfo)
    }
    # find lambda and/ or Awght
    if (!findAwght) {
      obj <- LKrigFindLambda(
        x = x,
        y = y,
        weights = weights,
        X = X,
        U = U,
        Z = Z,
        LKinfo = LKinfo,
        tol = tol,
        verbose = verbose
      )
      LKinfo <- LKinfoUpdate(LKinfo, lambda = obj$lambda.MLE)
    }
    else{
      obj <- LKrigFindLambdaAwght_hacked(
        x = x,
        y = y,
        X = X,
        weights = weights,
        U = U,
        Z = Z,
        LKinfo = LKinfo,
        verbose = verbose
      )
      LKinfo <- LKinfoUpdate(LKinfo,
                             lambda = obj$lambda.MLE,
                             a.wght = obj$a.wght.MLE)
    }
    # final call to LKrig to get all the summary statistics
    obj2 <- c(LKrig(
      x,
      y,
      weights = weights,
      Z = Z,
      X = X,
      U = U,
      LKinfo = LKinfo
    ),
    list(MLE = obj))
    class(obj2) <- c("LatticeKrig", "LKrig")
    obj2$call <- match.call()
    obj2$findAwght <- findAwght
    return(obj2)
  }
  
  
  
  LKrigFindLambdaAwght_hacked <- function(x,
                                          y,
                                          ...,
                                          LKinfo,
                                          use.cholesky = NULL,
                                          lowerBoundLogLambda =
                                            -16,
                                          upperBoundLogLambda = 4,
                                          lowerBoundOmega = -3,
                                          upperBoundOmega =  .75,
                                          factr = 1e7,
                                          pgtol = 1e-1,
                                          maxit = 15,
                                          verbose = FALSE) {
    #require(stats)
    
    # For rectangle omega = log(kappa) = log(sqrt(Awght-4))
    # but will change with other models.
    # Parts of the LKrig call that will be fixed.  (except updates to LKinfo)
    if (any(attr(LKinfo$a.wght, "isotropic")) == FALSE) {
      stop(
        paste(
          attr(LKinfo$a.wght, "isotropic"),
          "findAwght only setup to estimate a single a.wght
                 parameter in the model."
        )
      )
    }
    LKrigArgs <- c(list(x = x, y = y),
                   list(...),
                   list(
                     LKinfo = LKinfo,
                     NtrA = 0,
                     getVarNames = FALSE
                   ))
    
    if (verbose) {
      cat(
        "LKrigFindLambdaAwght: Set of LKrigArgs before first call:",
        names(LKrigArgs),
        fill = TRUE
      )
    }
    # Set up initial values of Awght and omega
    Awght.init <- as.numeric(LKrigArgs$LKinfo$a.wght[1])
    omega.start <- Awght2Omega(Awght.init, LKinfo)
    
    # Set up initial values of lambda and log lambda
    lambda.start <- LKrigArgs$LKinfo$lambda
    
    if (is.na(lambda.start)) {
      llambda.start <- -1
    }
    else{
      llambda.start <- log(lambda.start)
    }
    #
    if ((llambda.start < lowerBoundLogLambda) ||
        (llambda.start > upperBoundLogLambda) ||
        (is.na(llambda.start))) {
      stop("Given lambda value is out of bounds.")
    }
    #
    if (verbose) {
      cat(
        "LKrigFindLambdaAwght: llambda.start:",
        llambda.start,
        "a.wght.start:",
        Awght.init,
        fill = TRUE
      )
    }
    a.wghtTemp <- omega2Awght(omega.start, LKinfo)
    lambdaTemp <- exp(llambda.start)
    # initial call to likelihood and also to get symbolic decomposition of
    # the "M" matrix k-- sparsity pattern does not change as lambda, awght are varied.
    LKrigObject <- do.call("LKrig_hacked", c(
      LKrigArgs,
      list(
        use.cholesky = use.cholesky,
        return.cholesky = TRUE,
        return.wXandwU = TRUE,
        lambda = lambdaTemp,
        a.wght = a.wghtTemp,
        verbose = FALSE
      )
    ))
    # Update the LKrigArgs with cholesky decomp and wU
    LKrigArgs$use.cholesky <- LKrigObject$Mc
    LKrigArgs$wU <- LKrigObject$wUb
    #
    capture.evaluations <-  rbind(
      c(
        lambdaTemp,
        llambda.start,
        a.wghtTemp,
        omega.start,
        LKrigObject$sigma2.MLE.FULL,
        LKrigObject$tau.MLE.FULL,
        LKrigObject$lnProfileLike.FULL
      )
    )
    if (verbose) {
      cat("Capture.evaluations first call", fill = TRUE)
      cat(
        "lambda",
        "log lambda",
        "a.wght",
        "omega",
        "sigma2MLE",
        "tauMLE",
        "logProfileLike",
        fill = TRUE
      )
      cat(capture.evaluations, fill = TRUE)
    }
    
    #####  optimze likelihood over log lambda  and over omega =  log( a.wght -4)/2
    capture.env <- environment()
    # last two arguments are specific to this objectinve function
    result <- try(optim(
      c(llambda.start, omega.start),
      LambdaAwghtObjectiveFunction_hacked,
      lower = c(lowerBoundLogLambda, lowerBoundOmega),
      upper = c(upperBoundLogLambda, upperBoundOmega),
      method = "L-BFGS-B",
      #                      method="BFGS",
      control = list(
        fnscale = -1,
        factr = factr,
        pgtol = pgtol,
        maxit = maxit,
        ndeps = c(.05, .05)
      ),
      
      LKrigArgs = LKrigArgs,
      capture.env = capture.env,
      verbose = verbose
    ))
    if (verbose) {
      cat("Results from optimize:", fill = TRUE)
      print(result)
    }
    evalSummary <- !(class(result) == "try-error")
    llambda.MLE <- result$par[1]
    lambda.MLE <- exp(llambda.MLE)
    omega.MLE <- result$par[2]
    a.wght.MLE <- omega2Awght(omega.MLE, LKrigArgs$LKinfo)
    LKrigArgs$NtrA <- 20
    
    LKrigObject <- do.call("LKrig_hacked", c(LKrigArgs, list(
      lambda = lambda.MLE, a.wght = a.wght.MLE
    )))
    
    ###### end optimze block
    # save summary results from this set of parameters.
    # Output to be saved
    out <- rep(NA, 11)
    names(out) <-  c(
      "EffDf",
      "lnProfLike",
      "GCV",
      "tau.MLE",
      "sigma2.MLE",
      "lambda.MLE",
      "a.wght.MLE",
      "lnLike",
      "functionEval",
      "gradientEval",
      "totalEval"
    )
    out[1] <- LKrigObject$trA.est
    out[2] <- LKrigObject$lnProfileLike.FULL
    out[3] <- LKrigObject$GCV
    out[4] <- LKrigObject$tau.MLE.FULL
    out[5] <- LKrigObject$sigma2.MLE.FULL
    out[6] <- lambda.MLE
    out[7] <- a.wght.MLE
    out[8] <- LKrigObject$lnLike.FULL
    out[9] <- result$counts[1]
    out[10] <- result$counts[2]
    out[11] <- nrow(capture.evaluations)
    
    # Name columns  of likelihood eval.
    dimnames(capture.evaluations) <- list(
      NULL,
      c(
        "lambda",
        "logLambda",
        "a.wght",
        "omega",
        "sigma2.MLE",
        "tau.MLE",
        "lnProfileLike.FULL"
      )
    )
    return(
      list(
        summary = out,
        LKinfo = LKrigObject$LKinfo,
        llambda.start = llambda.start,
        Awght.start = omega2Awght(omega.start, LKrigArgs$LKinfo),
        lambda.MLE = lambda.MLE,
        a.wght.MLE = a.wght.MLE,
        omega.MLE = omega.MLE,
        llambda.MLE = llambda.MLE,
        lnLike.eval = capture.evaluations,
        call = match.call()
      )
    )
  }
  
  # Define the objective function
  LambdaAwghtObjectiveFunction_hacked <- function(PARS, LKrigArgs, capture.env, verbose =
                                                    FALSE) {
    lambdaTemp <- exp(PARS[1])
    a.wghtTemp <-  omega2Awght(PARS[2], LKrigArgs$LKinfo)
    
    hold <- do.call("LKrig_hacked", c(LKrigArgs, list(
      lambda = lambdaTemp, a.wght = a.wghtTemp
    )))[c("sigma2.MLE.FULL", "tau.MLE.FULL", "lnProfileLike.FULL")]
    rowForCapture <- c(
      lambdaTemp,
      PARS[1],
      a.wghtTemp,
      PARS[2],
      hold$sigma2.MLE.FULL,
      hold$tau.MLE.FULL,
      hold$lnProfileLike.FULL
    )
    if (verbose) {
      cat(rowForCapture, fill = TRUE)
    }
    lnProfileLike.FULL <- hold$lnProfileLike.FULL
    temp.eval <- get("capture.evaluations", envir = capture.env)
    assign("capture.evaluations",
           rbind(temp.eval, rowForCapture),
           envir = capture.env)
    return(lnProfileLike.FULL)
  }
  
  
  
  LKrig_hacked <- function(x,
                           y,
                           weights = NULL,
                           Z = NULL,
                           LKinfo = NULL,
                           iseed = NA,
                           NtrA = 20,
                           use.cholesky = NULL,
                           return.cholesky = TRUE,
                           X = NULL,
                           U = NULL,
                           wX = NULL,
                           wU = NULL,
                           return.wXandwU = TRUE,
                           ...,
                           getVarNames = TRUE,
                           verbose = FALSE) {
    # if LKinfo is missing create it from passed arguments
    # if it is passed update this object with the ... arguments
    # for example a new lambda value can be passed in this part.
    # LKinfo is a list that describes the LatticeKrig model.
    #   the list(...) only pertains to LKinfo arguments
    if (is.null(LKinfo)) {
      LKinfo <- do.call("LKrigSetup", c(list(x = x), list(...), list(verbose = verbose)))
    }
    else{
      LKinfo <- do.call("LKinfoUpdate", c(list(LKinfo = LKinfo), list(...)))
    }
    
    if (verbose) {
      cat(" ", fill = TRUE)
      cat("LKrig: updated LKinfo object", fill = TRUE)
      print(LKinfo)
      
    }
    #
    # the variable names are used to generate column labels if
    # those are missing.
    # the getVarNames switch is important -- if
    # LKrig is called via the do.call function then strange things
    # happen with the substitute function and I understand is an unresolved
    # aspect of Rli
    #
    if (getVarNames) {
      xName <- as.character(substitute(x))
      ZName <- as.character(substitute(Z))
      UName <- as.character(substitute(U))
      xName <- tail(xName, 1)
      ZName <- tail(ZName, 1)
      UName <- tail(UName, 1)
      
      # just take last component
      
    }
    else{
      xName <- "xVar"
      ZName <- "ZVar"
      UName <- "UVar"
    }
    #
    # create the initial parts of LKrig object
    # this list is added to as the computation proceeds
    # using the device  object<- list( object, newStuff)
    # and the full object is only obtained at the end
    # NOTE default for weights are just 1's and filled in by
    # the next call
    
    object <- createLKrigObject(
      x,
      y,
      weights = weights,
      Z = Z,
      X = X,
      U = U,
      LKinfo = LKinfo,
      xName = xName,
      ZName = ZName,
      UName = UName,
      verbose = verbose
    )
    
    nObs <-  nrow(object$y)
    nReps <- ncol(object$y)
    # for readablity make a local copy of LKinfo
    # but don't change it in this function!
    LKinfo <- object$LKinfo
    # Begin computations ....
    # weighted observation vector
    wy <- sqrt(object$weights) * object$y
    
    # create matrix for fixed part of model
    # Spatial drift matrix -- default is assumed to be linear in coordinates (m=2)
    # and includes possible covariate(s) -- the Z matrix.
    # the choice of fixed part of the model is controlled in LKinfo
    # (see also LKrigSetup)
    if (is.null(wU)) {
      wU <- LKrigMakewU(object, verbose = verbose)
    }
    
    # some column indices to keep track of fixed part of the model
    # NOTE nZ <= nt because Z is a subset of U
    object$nt <- ifelse(is.null(ncol(wU)), 0, ncol(wU))
    # create matrix for random part of model (basis functions)
    #  wX is the matrix of sum( N1*N2) basis function (columns) evaluated at the N locations (rows)
    # and multiplied by square root of diagonal weight matrix
    # this can be a large matrix if not encoded in sparse format.
    if (is.null(wX)) {
      timewX <- system.time(wX <- LKrigMakewX(object, verbose = verbose))
    }
    else{
      timewX <- rep(0, 5)
    }
    
    #   Precision matrix of the lattice process
    #   inverse of Q is proportional to the covariance matrix of the Markov Random Field
    timeQ <- system.time(Q <- nonstat_Q)
    # THIS IS HONESTLY THE ONLY BIG CHANGE I THINK
    # Q <- LKrig.precision(LKinfo, verbose=verbose)
    
    if (verbose) {
      print(LKinfo$a.wght[[1]] - 4)
    }
    # print(diag(Q))
    diag(Q) <- diag(Q) + (kappa2_weights * (LKinfo$a.wght[[1]] -
                                              4))
    # print(LKinfo$a.wght[[1]]-4) # for debugging
    
    if (LKinfo$dense) {
      if (!is.null(use.cholesky)) {
        stop("Can not update (use.cholesky) with dense matrix option")
      }
      Q <- spam2full(Q)
      wX <- spam2full(wX)
    }
    
    if (verbose & (!LKinfo$dense)) {
      cat("LKrig: Nonzero entries in Q:", length(Q@entries), fill = TRUE)
    }
    
    # G is the regularized (ridge) regression matrix that is
    # the key to the entire algorithm:
    # timeM<- system.time(
    timeM <- system.time(G <- t(wX) %*% wX + LKinfo$lambda * (Q))
    if (verbose) {
      if (!LKinfo$dense) {
        cat("LKrig: Nonzero entries in M:", length(G@entries), fill = TRUE)
      }
      else{
        cat("Dense matrix methods used", fill = TRUE)
      }
    }
    #  Find Cholesky square root of M
    #  This is where the heavy lifting happens!  M is in sparse, spam format so
    #  by the overloading this is a sparse cholesky decomposition.
    #  if this function has been coded efficiently this step should dominate
    #  all other computations ...
    #  If a previous sparse cholesky decoposition is passed then the
    #  pattern of sparseness is used for the decoposition.
    #  This can speed the computation as the symbolic decomposition part of the
    #  sparse Cholesky is a nontrivial step. The condition is that
    #  the current 'M' matrix  has the same sparse pattern as that
    #  which  will result in the same sparse pattern of factorization
    #  as when chol was applied to the case passed in 'use.cholesky'
    if (is.null(use.cholesky)) {
      timeChol <- system.time(GCholesky <- chol(G, memory = LKinfo$choleskyMemory))
    } else {
      timeChol <- system.time(GCholesky <- update.spam.chol.NgPeyton(use.cholesky, G))
    }
    
    if (!LKinfo$dense) {
      nonzero.entries <- length(GCholesky@entries)
    }
    else{
      nonzero.entries <- NA
    }
    
    
    if (verbose) {
      cat("LKrig: nonzero entries of GCholesky:",
          nonzero.entries,
          fill = TRUE)
    }
    # use GCholesky to find coefficients of estimate
    # Note that this functions also finds an important piece of the likelihood (quad.form)
    timeCoef <- system.time(
      out1 <- LKrig.coef(
        GCholesky,
        wX,
        wU,
        wy,
        LKinfo$lambda,
        collapseFixedEffect = LKinfo$collapseFixedEffect,
        verbose = verbose
      )
    )
    
    # Note collapseFixedEffect added as component here in the return
    # finding coefficients
    # fill in names of the fixed coefficients
    # fill in names of fixed model coefficients
    
    rownames(out1$d.coef) <- colnames(wU)
    #
    object <- c(object, out1)
    if (verbose) {
      cat("fixed model coefficients", fill = TRUE)
      cat(object$d.coef, fill = TRUE)
    }
    
    # compute predicted values  and residuals
    wfitted.values <- (wX %*% out1$c.coef)
    if (!is.null(wU)) {
      wfitted.values.fixed <- (wU %*% out1$d.coef)
      wfitted.values <- wfitted.values.fixed + wfitted.values
    }
    # X and U actully include the weights so need to divide these
    # out to get fitted values
    object$fitted.values <- wfitted.values / sqrt(object$weights)
    # For reference: fitted.values <- predict.LKrig(object, x, Znew = object$Z)
    # but at this point it is less efficient because X will be recomputed.
    object$residuals <- object$y - object$fitted.values
    # find likelihood
    timeLike <- system.time(
      #	out2 <- LKrig.lnPlikeOLD(GCholesky, Q, wy,
      #	             object$residuals, object$weights,
      #	             LKinfo)
      out2 <- LKrig.lnPlike(
        GCholesky,
        Q,
        object$quad.form,
        nObs,
        nReps,
        object$weights,
        LKinfo
      )
    )
    if (verbose) {
      cat("Likelihood/MLE list:", fill = TRUE)
      print(out2)
    }
    object <- c(object, out2)
    
    # estimate trace of hat matrix (effective degrees of freedom)
    # by Monte Carlo if NtrA greater than zero
    timeTrA <- system.time(if (NtrA > 0) {
      out3 <- LKrig.traceA(GCholesky,
                           wX,
                           wU,
                           LKinfo$lambda,
                           object$weights,
                           NtrA,
                           iseed = iseed)
      # find GCV using this trace estimate
      n <- length(object$weights)
      out3$GCV = (sum(object$weights * (object$residuals)^2) /
                    n) / (1 - out3$trA.est / n)^2
    } else {
      out3 <- list(trA.est = NA,
                   trA.SE = NA,
                   GCV = NA)
    })
    object <- c(object, out3)
    
    # create the table of times for individual function calls
    timingTable <- rbind(timewX, timeQ, timeM, timeChol, timeCoef, timeLike, timeTrA)
    timingTable <- timingTable[, 1:3]
    timingTable <- rbind(timingTable, colSums(timingTable))
    
    # last of required arguments to LKrig object
    object <- c(
      object,
      list(
        lambda.fixed = LKinfo$lambda,
        nonzero.entries = nonzero.entries,
        call = match.call(),
        timingLKrig = timingTable
      )
    )
    # finally add in some large matrices that can be reused if only
    # lambda is varied
    # NOTE to keep the code with the same components the name Mc is kept although
    # this is actually the cholesky decomposition of the G matrix as in the LKrig
    # article.
    if (return.cholesky) {
      object$Mc <- GCholesky
    }
    if (return.wXandwU) {
      object$wX <- wX
      object$wU <- wU
    }
    # refresh the class
    class(object) <- "LKrig"
    # all done!
    return(object)
  }
  
  
  
  ## import dataset ####
  station_time <- unique(stations_df$time)
  station_time <- station_time[order(station_time)]
  stations_df_day <- stations_df[stations_df$time == station_time[time_cycle], ]
  stations_df_day <- stations_df_day[!is.na(stations_df_day$EEA_NO2), ]
  stations_df_day[, "ssr"] <- stations_df_day[, "ssr"] / 10^6
  # skip days with less than 450 stations
  covars <- c("rh",
              "ssr",
              "t2m",
              "windspeed",
              "sl_blh",
              "EM_NO2",
              "DEM",
              "lag_cams_no2")
  if(ntry=="try7"){
    covars <- c("EM_NO2","DEM")
  }
  set.seed(1)
  N_data <- nrow(stations_df_day)
  ## start computation ####
  if (N_data >= 250) { #try3: from 450 to 250
    # & !paste0("data/output/10-fold/y_test_df_", time_cycle, ".rds") %in% list.files("data/output/10-fold")
    # print(paste0("Time cycle ", time_cycle, " already computed, skipping..."))
    if (surfaces) {
      ### 1a. Training data ####
      train_df <- stations_df_day
      s_train <- train_df[, c("Longitude", "Latitude")]
      y_train <- train_df$EEA_NO2
      Z_train <- as.matrix(train_df[, covars])
      
      #prediction grid
      pred_df[, "ssr"] <- pred_df[, "ssr"] / 10^6
      pred_time <- unique(pred_df$time)
      pred_time <- pred_time[order(pred_time)]
      pred_df_day <- pred_df[pred_df$time == pred_time[time_cycle], ]
      s_pred <- pred_df_day[, c("Longitude", "Latitude")]
      Z_pred <- as.matrix(pred_df_day[, covars])
      # Z_pred <- scale(Z_pred)
      pred_sp_day <- pred_df_day
      coordinates(pred_sp_day) <- c("Longitude", "Latitude")
      gridded(pred_sp_day) <- TRUE
      
      # 2. Models ####
      ## 2a. LM - Linear Model ####
      
      # training the model
      tic("time_lm")
      lm <- lm(
        EEA_NO2 ~ Longitude + Latitude + rh + ssr + t2m +
          windspeed + sl_blh + EM_NO2 + DEM + lag_cams_no2,
        data = train_df
      )
      time <- toc()$callback_msg
      
      # model interpretation
      # coefficients
      coef_ls <- as.data.frame(summary(lm)$coefficients)
      names(coef_ls) <- paste0("lm_", names(coef_ls))

      # surfaces
      pred_sp_day@data$lm <- predict(lm, newdata = pred_df_day)

      ## 2b. KED ####
      train_sp <- train_df
      coordinates(train_sp) <- c("Longitude","Latitude")
      # fit a variogram
      variogram <- variogram(EEA_NO2 ~ Longitude + Latitude + rh + ssr + t2m +
          windspeed + sl_blh + EM_NO2 + DEM + lag_cams_no2, data = train_sp)
      # plot(variogram)
      # # fit a variogram model
      variogram_model <- tryCatch(
        {
          fit.variogram(variogram, model = vgm("Exp"))
        },
        error = function(e) {
          print(paste("Errore:", e$message))
          return(NULL)
        }
      )

      if (is.null(variogram_model) | variogram_model$range[2]<0) {
      pred_sp_day@data$ked <- NA
      pred_sp_day@data$ked_sd <- NA
      } else {
      krige0 <- krige(formula = EEA_NO2 ~ Longitude + Latitude + rh + ssr + t2m +
            windspeed + sl_blh + EM_NO2 + DEM + lag_cams_no2,locations = train_sp,
                newdata=pred_sp_day,model=variogram_model)
              pred_sp_day@data$ked <- krige0$var1.pred
              # pred_sp_day@data$ked_sp <- pred_sp_day@data$ked - pred_sp_day@data$lm # qualcosa non va.. manca GLS?
              # pred_sp_day@data$ked_sd <- sqrt(krige0$var1.var)
              # g <- gstat(
              #             formula = EEA_NO2 ~ Longitude + Latitude + rh + ssr + t2m +
              #               windspeed + sl_blh + EM_NO2 + DEM + lag_cams_no2,
              #             data = train_sp,
              #             model = variogram_model
              #           )
        #       beta_gls <- predict(g, newdata = train_sp, BLUE = TRUE) # questo è il trend?
        # attr(beta_gls, "beta")
      }
      
      ## 2d. LK_stat - Lattice stationary ####
      ## 1 level of resolution, find Awght TRUE
      
      # training the model
      ##  remember -> new CTM grid:
      ##  drop latitude 35 - 35.1
      ##  and longitude 6 - 18.9
      rows_ctm <- length(seq(6.1, 18.8, by = .1))
      gridList_ctm <- list(x = seq(6.1, 18.8, by = .1),
                           y = seq(35.2, 47.9, by =
                                     .1))
      sGrid_ctm <- make.surface.grid(gridList_ctm)
      LKinfo_stat <- LKrigSetup(
        sGrid_ctm,
        a.wght = 4.01,
        # the rest of these you change directly in function calls
        NC = 128,
        nlevel = 1,
        normalize = TRUE,
        NC.buffer = 0
      )
      tic("time_LK_stat")
      LK_stat <- LatticeKrig(
        s_train,
        y_train,
        Z = Z_train,
        LKinfo = LKinfo_stat,
        findAwght = TRUE
      )
      time <- c(time, toc()$callback_msg)
      
      # coefficients
      # large scale
      coef_ls <- cbind(coef_ls, as.data.frame(LK_stat$d.coef))
      names(coef_ls)[ncol(coef_ls)] <- "LK_stat_Estimate"

      # spatial correlation
      sp_param <- data.frame(kappa2_stat = LK_stat$LKinfo$a.wght[[1]] - 4) # this is the kappa2 parameter that controls the spatial correlation in the model. The larger it is, the less spatial correlation there is.
    
      # nugget
      err_param <- data.frame(tau_stat = LK_stat$tau.MLE.FULL) # this is the nugget parameter that controls the amount of measurement error in the model. The larger it is, the more measurement error there is.
      err_param$sigma2_stat <- LK_stat$sigma2.MLE.FULL # this is the sigma2 parameter that controls the amount of spatial variability in the model. The larger it is, the more spatial variability there is.
      
      # add KED - try3
      if (!is.null(variogram_model)) {
      sp_param <- cbind(data.frame(range_ked=variogram_model$range[2],
                  sill_ked = variogram_model$psill[2]),sp_param)
      err_param <- cbind(data.frame(nugget_ked=variogram_model$psill[1]),err_param)
      } else {
        sp_param <- cbind(data.frame(range_ked=NA, nugget_ked=NA, sill_ked = NA),sp_param)
        err_param <- cbind(data.frame(nugget_ked=NA),err_param)
      }

      # surfaces
      if (!only_params){
      pred_sp_day@data$LK_stat <- predict(LK_stat, s_pred, Z = Z_pred)[, 1]
      pred_sp_day@data$LK_stat_ls <- predict(LK_stat,
                                             xnew = s_pred,
                                             Z = Z_pred,
                                             just.fixed = TRUE)[, 1]
      pred_sp_day@data$LK_stat_sp <- pred_sp_day@data$LK_stat - pred_sp_day@data$LK_stat_ls
      }

      # spplot(pred_sp_day, c("lm", "LK_stat_ls", "LK_stat_sp"))
      # spplot(pred_sp_day, c("LK_stat_sp"))
      
      ## 2e. LK_stun - Lattice Kriging + STUN parameters ####
      
      # import STUN parameters
      # params <- params_all[, 128:1, , time_cycle] # which params we want?
      # names_param <- c("surr30", "surr1", "rscale30", "rscale1")

      # for (i in 1) # loop through the 4 different STUN parameter sets, 1: surr30, 2: surr1, 3: rscale30, 4: rscale1
      # {
      #   if (i == 1){
      #     params <- params_all_surr30[, 128:1, , time_cycle] # which params we want?
      #   } else if (i == 2){
      #     params <- params_all_surr1[, 128:1, , time_cycle] # which params we want?
      #   } else if (i == 3){
      #     params <- params_all_rscale30[, 128:1, , time_cycle] # which params we want?
      #   } else {
      #     params <- params_all_rscale1[, 128:1, , time_cycle] # which params we want?
      #   }
      params_all <- params_all_surr30
      params <- params_all[, 128:1, , time_cycle]
      kappa2 <- exp(params[, , 1])
      awght <- kappa2 + 4
      theta <- params[, , 2] + pi / 2
      rho   <- params[, , 3]
      rhox <- sqrt(rho)
      rhoy <- 1 / rhox
      H11 <- (rhox^2 * (cos(theta))^2) + (rhoy^2 * (sin(theta))^2)
      H12 <- (rhoy^2 - rhox^2) * (sin(theta) * cos(theta))
      H21 <- H12
      H22 <- (rhox^2 * (sin(theta))^2) + (rhoy^2 * (cos(theta))^2)
      stencil_tensor <- array(NA, c(rows_ctm, rows_ctm, 9))
      stencil_tensor[, , 1] <- 0.5 * H12
      stencil_tensor[, , 2] <- -H22
      stencil_tensor[, , 3] <- -0.5 * H12
      stencil_tensor[, , 4] <- -H11
      stencil_tensor[, , 5] <- kappa2 + 2 * H11 + 2 * H22
      stencil_tensor[, , 6] <- -H11
      stencil_tensor[, , 7] <- -0.5 * H12
      stencil_tensor[, , 8] <- -H22
      stencil_tensor[, , 9] <- 0.5 * H12
      
      awght_obj <- list(x = gridList_ctm$x, y = gridList_ctm$y, z = stencil_tensor) #gridList_ctm from LK_stun section
      class(awght_obj) <- "multivariateSurfaceGrid"
      
      # training the model
      LKinfo_stun <- LKrigSetup(
        sGrid_ctm,
        a.wghtObject = awght_obj,
        NC = 128,
        nlevel = 1,
        normalize = T,
        NC.buffer = 0
      )
      nonstat_Q <- LKrig.precision(LKinfo_stun)
      # tic("time_LK_stun")
      # LK_stun <- LatticeKrig(
      #   x = s_train,
      #   y = y_train,
      #   Z = Z_train,
      #   LKinfo = LKinfo_stun
      # )
      # time <- c(time, toc()$callback_msg)
      # # surfaces
      # pred_sp_day@data$LK_stun <- predict(LK_stun, s_pred, Z = Z_pred)[, 1]
      # pred_sp_day@data$LK_stun_ls <- predict(LK_stun,
      #                                        xnew = s_pred,
      #                                        Z = Z_pred,
      #                                        just.fixed = TRUE)[, 1]
      # pred_sp_day@data$LK_stun_sp <- pred_sp_day@data$LK_stun - pred_sp_day@data$LK_stun_ls
      
      
      ## 2f. LK_stun_adj - Adjusted Lattice Krig + STUN ####
      
      #training the model
      # kappa2_weights <- 1
      
      tic("time_LK_stun_adj")
      LKinfo_stun_adj <- est_additional_kappa2(s_train, y_train, Z_train, LKinfo_stun, kappa2_weights = kappa2_weights)
      LK_stun_adj <- LatticeKrig(
        x = s_train,
        y = y_train,
        Z = Z_train,
        LKinfo = LKinfo_stun_adj
      )
      time <- c(time, toc()$callback_msg)
      
      # coefficients
      # large scale
      coef_ls <- cbind(coef_ls, as.data.frame(LK_stun_adj$d.coef))
      names(coef_ls)[ncol(coef_ls)] <- paste("LK_stun_adj_Estimate")
      
      # spatial correlation
      kappa2_adj <- LK_stun_adj$LKinfo$a.wght[[1]][,5] + (2 * LK_stun_adj$LKinfo$a.wght[[1]][,2]) + (2 * LK_stun_adj$LKinfo$a.wght[[1]][,4])
      kappa2_stun <- LKinfo_stun$a.wght[[1]][,5] + (2 * LKinfo_stun$a.wght[[1]][,2]) + (2 * LKinfo_stun$a.wght[[1]][,4])
      kappa2_eea <- kappa2_adj - kappa2_stun
      sp_param <- cbind(sp_param,data.frame(kappa2_adj=kappa2_adj,kappa2_stun=kappa2_stun,kappa2_eea=kappa2_eea))  # this is the kappa2 parameter that controls the spatial correlation in the model. The larger it is, the less spatial correlation
      sp_param <- cbind(sp_param, data.frame(kappa2_weights = kappa2_weights)) # this is the kappa2 parameter that controls the spatial correlation in the model. The larger it is, the less spatial correlation
       
      # nugget effect
      err_param <- cbind(err_param, data.frame(tau_stun_adj = LK_stun_adj$tau.MLE.FULL)) # this is the nugget parameter that controls the amount of measurement error in the model. The larger it is, the more measurement error there is.
      err_param <- cbind(err_param, data.frame(sigma2_stun_adj = LK_stun_adj$sigma2.MLE.FULL)) # this is the sigma2 parameter that controls the amount of spatial variability in the model. The larger it is, the more spatial variability there is.
      
      # surfaces
      if (!only_params){
      pred_sp_day@data$LK_stun_adj <- predict(LK_stun_adj, s_pred, Z = Z_pred)[, 1]
      pred_sp_day@data$LK_stun_adj_ls <- predict(LK_stun_adj,
                                                 xnew = s_pred,
                                                 Z = Z_pred,
                                                 just.fixed = TRUE)[, 1]
      pred_sp_day@data$LK_stun_adj_sp <- pred_sp_day@data$LK_stun_adj - pred_sp_day@data$LK_stun_adj_ls
      }
      # 3. exporting ####
      if (!only_params){
      save(pred_sp_day,
           file = file.path(map_folder,
            paste0("prediction_map_", time_cycle, ".rda"))) #10 MB
      }
      coef_ls$time <- station_time[time_cycle]
      sp_param$time <- station_time[time_cycle]
      err_param$time <- station_time[time_cycle]
      save(coef_ls,sp_param,err_param,file = file.path(coef_folder,
        paste0("coef_", time_cycle, ".rda")))
      # ADD TIME !!
      
      # spplot(pred_sp_day,"LK_stun_adj_sp")
      # spplot(pred_sp_day,c("ked_sp","LK_def_sp","LK_stat_sp","LK_stun_sp","LK_stun_adj_sp"))
      # pdf(paste0("data/output/plot/results_", time_cycle, ".pdf"))
      # print(spplot(pred_sp_day, c("LK_stat_sp", "LK_stun_adj_sp")))
      # print(spplot(pred_sp_day, c("LK_stat_ls", "LK_stun_adj_ls")))
      # print(spplot(pred_sp_day, c("lm", "LK_stat", "LK_stun_adj")))
      # dev.off()
      # try(dev.off())
      # coef_ls$time <- station_time[time_cycle]
      # sp_param$time <- station_time[time_cycle]
      # save(coef_ls, sp_param, file = paste0("data/output/coef/coef_sp_param_", time_cycle, ".rda"))
    }
    
    if (CV_scheme) {
      kfold <- list()
      
      # big-regions
      if (CV_type == "big-regions") {
        step_kfold <- c(1, N_data * seq(0.1, 1, 0.1))
        cv_idx <- sample(N_data, N_data, replace = F)
        for (nk in 1:10) {
          start_fold <- ceiling(step_kfold[nk])
          end_fold <- floor(step_kfold[nk + 1])
          kfold[[nk]] <- cv_idx[start_fold:end_fold]
        }
      }
      
      # 10-fold
      if (CV_type == "10-fold") {
        step_kfold <- c(1, N_data * seq(0.1, 1, 0.1))
        cv_idx <- sample(N_data, N_data, replace = F)
        for (nk in 1:10) {
          start_fold <- ceiling(step_kfold[nk])
          end_fold <- floor(step_kfold[nk + 1])
          kfold[[nk]] <- cv_idx[start_fold:end_fold]
        }
      }
      for (nk in 1:length(kfold)) {
        test_ind <- kfold[[nk]]
        ## 1a. Training data ####
        #  cams_test <- sample(unique(stations_df_day$CAMS_NO2), length(unique(stations_df_day$CAMS_NO2)) /
        #                       10)
        # test_ind <- which(stations_df_day$CAMS_NO2 %in% cams_test)
        # test_ind <- sample(1:nrow(stations_df_day), floor(nrow(stations_df_day) /
        #                                                     10))
        # stations_df_day[,c("Longitude","Latitude",covars)] <- scale(stations_df_day[,c("Longitude","Latitude",covars)])
        train_df <- stations_df_day[-test_ind, ]
        test_df  <- stations_df_day[test_ind, ]
        
        s_train <- train_df[, c("Longitude", "Latitude")]
        y_train <- train_df$EEA_NO2
        Z_train <- as.matrix(train_df[, covars])
        s_test <- test_df[, c("Longitude", "Latitude")]
        y_test <- test_df$EEA_NO2
        Z_test <- as.matrix(test_df[, covars])
        # Z_test <- scale(Z_test)
        
        ## 1b. Prediction grid ####
        # pred_df[, "ssr"] <- pred_df[, "ssr"] / 10^6
        #
        # pred_time <- unique(pred_df$time)
        # pred_time <- pred_time[order(pred_time)]
        # pred_df_day <- pred_df[pred_df$time == pred_time[time_cycle], ]
        #
        # s_pred <- pred_df_day[, c("Longitude", "Latitude")]
        # Z_pred <- as.matrix(pred_df_day[, covars])
        # # Z_pred <- scale(Z_pred)
        #
        # pred_sp_day <- pred_df_day
        # coordinates(pred_sp_day) <- c("Longitude", "Latitude")
        # gridded(pred_sp_day) <- TRUE
        
        # 2. Models ####
        names_test_df <- c("lwr", "fit", "upr")
        
        ## 2a. LM - Linear Model ####
        
        # training the model
        tic("time_lm")
        lm <- lm(
          EEA_NO2 ~ Longitude + Latitude + rh + ssr + t2m +
            windspeed + sl_blh + EM_NO2 + DEM + lag_cams_no2,
          data = train_df
        )
        time <- toc()$callback_msg
        
        # surfaces
        # pred_sp_day@data$lm <- predict(lm, newdata = pred_df_day)
        
        # test predictions
        y_test_lm <- predict(lm, newdata = test_df, interval = "prediction")
        y_test_lm <- y_test_lm[, c(2, 1, 3)]
        colnames(y_test_lm) <- paste0("lm_", colnames(y_test_lm))
        y_test_df <- as.data.frame(y_test_lm)
        
        ## 2b. KED - Kriging with External Drift (KED) ####
        
        # training the model
        # tic("time_ked")
        train_sp <- train_df
        coordinates(train_sp) <- c("Longitude","Latitude")
        test_sp <- test_df
        coordinates(test_sp) <- c("Longitude","Latitude")
        tic("time_KED")
        variogram <- variogram(EEA_NO2 ~ Longitude + Latitude + rh + ssr + t2m +
          windspeed + sl_blh + EM_NO2 + DEM + lag_cams_no2, data = train_sp)
      # plot(variogram)
      # # fit a variogram model
      variogram_model <- tryCatch(
        {
          fit.variogram(variogram, model = vgm("Exp"))
        },
        error = function(e) {
          print(paste("Errore:", e$message))
          return(NULL)
        }
      )


      if (is.null(variogram_model) | variogram_model$range[2]<0) {
      y_test_ked <- rep(NA,length(test_sp))
      se_ked <- rep(NA,length(test_sp))
      } else {
      krige0 <- krige(formula = EEA_NO2 ~ Longitude + Latitude + rh + ssr + t2m +
            windspeed + sl_blh + EM_NO2 + DEM + lag_cams_no2,locations = train_sp,
                newdata=test_sp,model=variogram_model)
        y_test_ked <- krige0$var1.pred
        se_ked <- sqrt(krige0$var1.var)
      }
      
        time <- c(time, toc()$callback_msg)
        y_test_ked <- cbind(
          y_test_ked - 1.96 * se_ked,
          y_test_ked,
          y_test_ked + 1.96 * se_ked
        )
        colnames(y_test_ked) <- names_test_df
        colnames(y_test_ked)<-paste0("ked_",colnames(y_test_ked))
        y_test_df <- cbind(y_test_df,y_test_ked)
        
        
        # ked <- spatialProcess(s_train, y_train, Z = Z_train)
        
        # print("error")
        # print(time_cycle)
        # time <- c(time,toc()$callback_msg)
        
        
        # surfaces
        # pred_sp_day@data$ked <- predict(ked, s_pred, Z = Z_pred)
        # pred_sp_day@data$ked_ls <- predict(ked, xnew=s_pred, Z=Z_pred, just.fixed = TRUE)
        # pred_sp_day@data$ked_sp <- pred_sp_day@data$ked - pred_sp_day@data$ked_ls
        # spplot(pred_sp_day,"ked_sp")
        
        # # test predictions
        # y_test_ked <- predict(ked, s_test, Z = Z_test)
        # se_ked <- predictSE(ked, s_test, Z = Z_test)
        # y_test_ked <- cbind(
        #   y_test_ked - 1.96 * sqrt(se_ked^2 + ked[["summary"]][["tau"]]^2),
        #   y_test_ked,
        #   y_test_ked + 1.96 * sqrt(se_ked^2 + ked[["summary"]][["tau"]]^2)
        # )
        # colnames(y_test_ked) <- names_test_df
        # colnames(y_test_ked)<-paste0("ked_",colnames(y_test_ked))
        # y_test_df <- cbind(y_test_df,y_test_ked)
        
        ## 2c. LK_def - Lattice Krig Default ####
        
        # training the model
        tic("time_LK_def")
        LK_def <- LatticeKrig(s_train, y_train, Z = Z_train)
        time <- c(time, toc()$callback_msg)
        
        # surfaces
        # pred_sp_day@data$LK_def <- predict(LK_def, s_pred, Z = Z_pred)
        # pred_sp_day@data$LK_def_ls <- predict(LK_def,
        #                                       xnew = s_pred,
        #                                       Z = Z_pred,
        #                                       just.fixed = TRUE)
        # pred_sp_day@data$LK_def_sp <- pred_sp_day@data$LK_def - pred_sp_day@data$LK_def_ls
        # spplot(pred_sp_day,"LK_def_sp")
        
        # test predictions
        y_test_LK_def <- predict(LK_def, s_test, Z = Z_test)
        # y_test_LK_def_ls <- predict(LK_def, s_test, Z = Z_test, just.fixed = TRUE)
        se_LK_def <- predictSE(LK_def, s_test, Z = Z_test)
        y_test_LK_def <- cbind(
          y_test_LK_def - 1.96 * sqrt(se_LK_def^2 + LK_def$tau.MLE^2),
          y_test_LK_def,
          y_test_LK_def + 1.96 * sqrt(se_LK_def^2 + LK_def$tau.MLE^2)
        )
        colnames(y_test_LK_def) <- names_test_df
        colnames(y_test_LK_def) <- paste0("LK_def_", colnames(y_test_LK_def))
        y_test_df <- cbind(y_test_df, y_test_LK_def)
        
        ## 2d. LK_stat - Lattice stationary ####
        ## 1 level of resolution, find Awght TRUE
        
        # training the model
        ##  remember -> new CTM grid:
        ##  drop latitude 35 - 35.1
        ##  and longitude 6 - 18.9
        rows_ctm <- length(seq(6.1, 18.8, by = .1))
        gridList_ctm <- list(x = seq(6.1, 18.8, by = .1),
                             y = seq(35.2, 47.9, by =
                                       .1))
        sGrid_ctm <- make.surface.grid(gridList_ctm)
        LKinfo_stat <- LKrigSetup(
          sGrid_ctm,
          a.wght = 4.01,
          # the rest of these you change directly in function calls
          NC = 128,
          nlevel = 1,
          normalize = TRUE,
          NC.buffer = 0
        )
        tic("time_LK_stat")
        LK_stat <- LatticeKrig(
          s_train,
          y_train, 
          Z = Z_train,
          LKinfo = LKinfo_stat,
          findAwght = TRUE
        )
        time <- c(time, toc()$callback_msg)
        
        # surfaces
        # pred_sp_day@data$LK_stat <- predict(LK_stat, s_pred, Z = Z_pred)[, 1]
        # pred_sp_day@data$LK_stat_ls <- predict(LK_stat,
        #                                        xnew = s_pred,
        #                                        Z = Z_pred,
        #                                        just.fixed = TRUE)[, 1]
        # pred_sp_day@data$LK_stat_sp <- pred_sp_day@data$LK_stat - pred_sp_day@data$LK_stat_ls
        # # spplot(pred_sp_day,"LK_stat_sp")
        
        # test predictions
        y_test_LK_stat <- predict(LK_stat, s_test, Z = Z_test)
        se_LK_stat <- predictSE(LK_stat, s_test, Z = Z_test)
        y_test_LK_stat <- cbind(
          y_test_LK_stat - 1.96 * sqrt(se_LK_stat^2 + LK_stat$tau.MLE^2),
          y_test_LK_stat,
          y_test_LK_stat + 1.96 * sqrt(se_LK_stat^2 + LK_stat$tau.MLE^2)
        )
        colnames(y_test_LK_stat) <- names_test_df
        colnames(y_test_LK_stat) <- paste0("LK_stat_", colnames(y_test_LK_stat))
        y_test_df <- cbind(y_test_df, y_test_LK_stat)
        
        ## 2e. LK_stun - Lattice Kriging + STUN parameters ####
        names_param <- c("surr30", "surr1", "rscale30", "rscale1")
        # import STUN parameters
        for (i in 1) # loop through the 4 different STUN parameter sets, 1: surr30, 2: surr1, 3: rscale30, 4: rscale1
        {
          if (i == 1){
            params <- params_all_surr30[, 128:1, , time_cycle] # which params we want?
          } else if (i == 2){
            params <- params_all_surr1[, 128:1, , time_cycle] # which params we want?
          } else if (i == 3){
            params <- params_all_rscale30[, 128:1, , time_cycle] # which params we want?
          } else {
            params <- params_all_rscale1[, 128:1, , time_cycle] # which params we want?
          }
        
        #params <- params_all[, 128:1, , time_cycle]
        kappa2 <- exp(params[, , 1])
        awght <- kappa2 + 4
        theta <- params[, , 2] + pi / 2
        rho   <- params[, , 3]
        rhox <- sqrt(rho)
        rhoy <- 1 / rhox
        H11 <- (rhox^2 * (cos(theta))^2) + (rhoy^2 * (sin(theta))^2)
        H12 <- (rhoy^2 - rhox^2) * (sin(theta) * cos(theta))
        H21 <- H12
        H22 <- (rhox^2 * (sin(theta))^2) + (rhoy^2 * (cos(theta))^2)
        stencil_tensor <- array(NA, c(rows_ctm, rows_ctm, 9))
        stencil_tensor[, , 1] <- 0.5 * H12
        stencil_tensor[, , 2] <- -H22
        stencil_tensor[, , 3] <- -0.5 * H12
        stencil_tensor[, , 4] <- -H11
        stencil_tensor[, , 5] <- kappa2 + 2 * H11 + 2 * H22
        stencil_tensor[, , 6] <- -H11
        stencil_tensor[, , 7] <- -0.5 * H12
        stencil_tensor[, , 8] <- -H22
        stencil_tensor[, , 9] <- 0.5 * H12
        
        awght_obj <- list(x = gridList_ctm$x,
                          y = gridList_ctm$y,
                          z = stencil_tensor) #gridList_ctm from LK_stun section
        class(awght_obj) <- "multivariateSurfaceGrid"
        
        # training the model
        LKinfo_stun <- LKrigSetup(
          sGrid_ctm,
          a.wghtObject = awght_obj,
          NC = 128,
          nlevel = 1,
          normalize = T,
          NC.buffer = 0
        )
        nonstat_Q <- LKrig.precision(LKinfo_stun)
        tic("time_LK_stun")
        LK_stun <- LatticeKrig(
          x = s_train,
          y = y_train,
          Z = Z_train,
          LKinfo = LKinfo_stun
        )
        time <- c(time, toc()$callback_msg)
        
        # surfaces
        # pred_sp_day@data$LK_stun <- predict(LK_stun, s_pred, Z = Z_pred)[, 1]
        # pred_sp_day@data$LK_stun_ls <- predict(LK_stun,
        #                                        xnew = s_pred,
        #                                        Z = Z_pred,
        #                                        just.fixed = TRUE)[, 1]
        # pred_sp_day@data$LK_stun_sp <- pred_sp_day@data$LK_stun - pred_sp_day@data$LK_stun_ls
        # spplot(pred_sp_day,"LK_stun_sp")
        
        # test predictions
        y_test_LK_stun <- predict(LK_stun, s_test, Z = Z_test)
        se_LK_stun <- predictSE(LK_stun, s_test, Z = Z_test)
        y_test_LK_stun <- cbind(
          y_test_LK_stun - 1.96 * sqrt(se_LK_stun^2 + LK_stun$tau.MLE^2),
          y_test_LK_stun,
          y_test_LK_stun + 1.96 * sqrt(se_LK_stun^2 + LK_stun$tau.MLE^2)
        )
        colnames(y_test_LK_stun) <- names_test_df
        colnames(y_test_LK_stun) <- paste0("LK_stun_", colnames(y_test_LK_stun))
        y_test_df <- cbind(y_test_df, y_test_LK_stun)
        
        ## 2f. LK_stun_adj - Adjusted Lattice Krig + STUN ####
        
        #training the model
        # kappa2_weights <- 1
        tic("time_LK_stun_adj")
        LKinfo_stun_adj <- est_additional_kappa2(s_train,
                                                 y_train,
                                                 Z_train,
                                                 LKinfo_stun,
                                                 kappa2_weights = kappa2_weights)
        LK_stun_adj <- LatticeKrig(
          x = s_train,
          y = y_train,
          Z = Z_train,
          LKinfo = LKinfo_stun_adj
        )
        time <- c(time, toc()$callback_msg)
        
        # surfaces
        # pred_sp_day@data$LK_stun_adj <- predict(LK_stun_adj, s_pred, Z = Z_pred)[, 1]
        # pred_sp_day@data$LK_stun_adj_ls <- predict(LK_stun_adj,
        #                                            xnew = s_pred,
        #                                            Z = Z_pred,
        #                                            just.fixed = TRUE)[, 1]
        # pred_sp_day@data$LK_stun_adj_sp <- pred_sp_day@data$LK_stun_adj - pred_sp_day@data$LK_stun_adj_ls
        # spplot(pred_sp_day,"LK_stun_adj_sp")
        # spplot(pred_sp_day,c("ked_sp","LK_def_sp","LK_stat_sp","LK_stun_sp","LK_stun_adj_sp"))
        
        # test predictions
        y_test_LK_stun_adj <- predict(LK_stun_adj, s_test, Z = Z_test)
        se_LK_stun_adj <- predictSE(LK_stun_adj, s_test, Z = Z_test)
        y_test_LK_stun_adj <- cbind(
          y_test_LK_stun_adj - 1.96 * sqrt(se_LK_stun_adj^2 + LK_stun_adj$tau.MLE^2),
          y_test_LK_stun_adj,
          y_test_LK_stun_adj + 1.96 * sqrt(se_LK_stun_adj^2 + LK_stun_adj$tau.MLE^2)
        )
        colnames(y_test_LK_stun_adj) <- names_test_df
        colnames(y_test_LK_stun_adj) <- paste0("LK_stun_adj_",colnames(y_test_LK_stun_adj))
        # colnames(y_test_LK_stun_adj) <- paste0("LK_stun_adj_", colnames(y_test_LK_stun_adj))
        y_test_df <- cbind(y_test_df, y_test_LK_stun_adj)
        }
        # ## 2g. FRK - Fixed Rank Kriging ####
        #
        # # training model
        # df_FRK <- as.data.frame(cbind(y_train, s_train))
        # names(df_FRK)[1] <- "AQ_EEA_NO2"
        # f_FRK <- as.formula(paste("AQ_EEA_NO2~", paste(
        #   c(
        #     "rh",
        #     "ssr",
        #     "t2m",
        #     "windspeed",
        #     "sl_blh",
        #     "EM_NO2",
        #     "lag_cams_no2"
        #   ),
        #   collapse = "+"
        # )))
        # sp_FRK <- df_FRK
        # coordinates(sp_FRK) <- c("Longitude", "Latitude")
        # tic("FRK")
        # FRK <- FRK(f_FRK, sp_FRK, BAUs = pred_sp_day)
        # time <- c(time,toc()$callback_msg)
        #
        # # surfaces
        # pred_grid_FRK <-
        #   predict(FRK, obs_fs = FALSE)
        # pred_sp_day@data$FRK <- pred_grid_FRK@data$mu
        # pred_sp_day@data$FRK_ls <- (cbind(1,Z_pred) %*% FRK@alphahat)[,1]
        # pred_sp_day@data$FRK_sp <- (FRK@S0 %*% FRK@mu_eta)[,1]
        # pred_sp_day@data$FRK_xi <- pred_sp_day@data$FRK - (pred_sp_day@data$FRK_ls + pred_sp_day@data$FRK_sp)
        # # spplot(pred_sp_day,c("ked_sp","LK_def_sp","LK_stat_sp","LK_stun_sp","LK_stun_adj_sp","FRK_sp"))
        # # spplot(pred_sp_day,c("ked_ls","LK_def_ls","LK_stat_ls","LK_stun_ls","LK_stun_adj_ls","FRK_ls"))
        # # spplot(pred_sp_day,c("lm","ked","LK_def","LK_stat","LK_stun","LK_stun_adj","FRK"))
        # # spplot(pred_sp_day,c("FRK_xi"))
        #
        # # test predictions
        # points_test <- s_test
        # coordinates(points_test) <- c("Longitude", "Latitude")
        # points_test <- SpatialPointsDataFrame(points_test, as.data.frame(y_test))
        # y_test_FRK <- cbind(points_test@data, over(points_test, pred_grid_FRK))[, c("mu", "sd")]
        # y_test_FRK <- cbind(
        #   y_test_FRK$mu - 1.96 * sqrt(y_test_FRK$sd^2 + FRK@sigma2fshat),
        #   y_test_FRK$mu,
        #   y_test_FRK$mu + 1.96 * sqrt(y_test_FRK$sd^2 + FRK@sigma2fshat)
        # )
        # colnames(y_test_FRK) <- names_test_df
        # colnames(y_test_FRK)<-paste0("FRK_",colnames(y_test_FRK))
        # y_test_df <- cbind(y_test_df,y_test_FRK)
        # y_test_df <- cbind(y_test,y_test_df)
        
        # 3. Exporting outputs ####
        
        # models
        # save(lm,ked,LK_def,LK_stat,LK_stun,LK_stun_adj,file=paste0("data/output/models_",time_cycle,".rda")) # 53 MB for each day
        
        # surfaces
        #  save(pred_sp_day,file = paste0("data/output/prediction_map_",time_cycle,".rda")) #10 MB
        y_test_df <- cbind(test_df[, c("AirQualityStation", "time")], y_test_df)
        # test predictions
        if (nk == 1) {
          all_y_test_df <- y_test_df
        } else{
          all_y_test_df <- rbind(all_y_test_df, y_test_df)
        }
        
        # time
        # saveRDS(time, file = paste0("data/output/time_", time_cycle, ".rds")) # 7.8 KB
        
        # pdf(paste0("plot/results", time_cycle, ".pdf"))
        # spplot(pred_sp_day,
        #        c("LK_def_sp", "LK_stat_sp", "LK_stun_sp", "LK_stun_adj_sp"))
        # spplot(pred_sp_day,
        #        c("LK_def_ls", "LK_stat_ls", "LK_stun_ls", "LK_stun_adj_ls"))
        # spplot(pred_sp_day,
        #        c("lm", "LK_def", "LK_stat", "LK_stun", "LK_stun_adj"))
        # dev.off()
        
      }
      # CV_folder
      saveRDS(all_y_test_df,
              file = file.path(CV_folder,paste0("y_test_df_", time_cycle, ".rds"))) # 7.8 KB
    }
  }
}


  