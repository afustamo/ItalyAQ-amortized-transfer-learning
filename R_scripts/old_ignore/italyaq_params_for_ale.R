# -----------------------------------------------------------------------
# INFO: 
# -----------------------------------------------------------------------

# I converted the data you sent me into 
# an h5 file called "IAQ_NO2_resid_2019.h5"
# I used this data to create these parameter datasets 
# with a naming convention: 

# -----------------------------------------------------------------------
# NAMING CONVENTION 
# -----------------------------------------------------------------------

# 1 rep means only that exact field was used
# 30 reps means that the 14 before and 15 after were used

# surround means the 14 before and 15 after were also used to normalize
# that field

# rscale means that the individual field was z score normalized 

# here we open up df_rscale_30rep.h5 
# (this is probably the one you guys should use anyway)

# -----------------------------------------------------------------------
# DATA EXTRACTION
# ----------------------------------------------------------------------- 

# here is how i install the rhdf5 package
# i also sent u a file on how to work with h5 called "h5_for_ale.R"
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("rhdf5")

#libraries 
library(rhdf5)
library(LatticeKrig)
library(maps)

# load in actual clim fields
AQ_file_path <- "data/IAQ_NO2_resid_2019.h5"

# load in actual clim fields (these are in data/ in the google drive)
AQ_file_path <- "C:/Users/anton/Desktop/Research/Italy_AQ_paper/Italy_AQ_AmortisedLatticeKrig/data/IAQ_NO2_resid_2019.h5"

df_AQ <- h5read(AQ_file_path, "fields")
lat <- rev(h5read(AQ_file_path, "lat"))
lon <- rev(h5read(AQ_file_path, "lon"))
time <- h5read(AQ_file_path, "time")


# load in i2i network clim outputs
file_path <- "data/df_rscale_30rep.h5"

# load in i2i network clim outputs (these are in results/ in the google drive)
file_path <- "C:/Users/anton/Desktop/Research/Italy_AQ_paper/Italy_AQ_AmortisedLatticeKrig/results/italy_aq_outputs/df_rscale_30rep.h5"
h5ls(file_path)

# parameters from global stun estimator
df_STUN <- h5read(file_path, "stun")

# parameters from local cnn estimator (we wont use this now)
df_CNN <- h5read(file_path, "cnn")

# params should have dim 128 x 128 x 3 x 100
# this is 100 days (like you sent me)
# 128 x 128 is the field size
# in dimension 3 the first field is log(kappa2)
# the second is theta, the third is rho 


# extract/make params for STUN
# we will only work with these in this file because the CNN ones are garbage
kappa2_s <- exp(df_STUN[,,1,])
awght_s <- kappa2_s + 4
# need to transform theta for LK 
theta_s <- df_STUN[,,2,] + pi/2
# for plotting intuitive angles
rho_s <- df_STUN[,,3,]


# also extract/make params for CNN
# these seem to be garbage though so i wouldn't use 
kappa2_c <- exp(df_CNN[,,1,])
awght_c <- kappa2_c + 4
# yes, this theta needs to be transformed differently 
theta_c <- ((-1.0) *  df_CNN[,,2,])
rho_c <- df_CNN[,,3,]


# shift everything to look nice on plots
df_AQ <- aperm(df_AQ, c(2,1,3))
df_AQ <- df_AQ[, dim(df_AQ)[2]:1, ]
kappa2_s <- kappa2_s[, dim(kappa2_s)[2]:1, ]
theta_s <- theta_s[, dim(theta_s)[2]:1, ]
rho_s <- rho_s[, dim(rho_s)[2]:1, ]


# for later coding into LK
rhox_s <- sqrt(rho_s)
rhoy_s <- 1/rhox_s


# -----------------------------------------------------------------------
# PLOTTING DATA AND PARAMS
# -----------------------------------------------------------------------


# plotting the first field and it's associated params
library(maps)
plot_params <- function(
    idx, 
    aq_fields, 
    kappa2,
    theta, 
    rho, 
    border_lwd = 1
){
  awght <- kappa2[,,idx] + 4
  # sanity plotting for input and output params
  par(mfrow = c(2,2), mar =   c(3.1, 4.1, 2.1, 2.1))
  imagePlot(x = lon, y = lat, 
            aq_fields[,,idx], main = "NO2 Residual Field", col = turbo(256)) #, 
  #horizontal = TRUE)
  map("world2", add = TRUE, 
      col = "grey0", lwd = border_lwd)
  
  imagePlot(x = lon, y = lat,
            awght, main = "awght", col = viridis(256))
  map("world2", add = TRUE, 
      col = "grey90", lwd = border_lwd)
  
  imagePlot(x = lon, y = lat,
            theta[,,idx], main = "Theta", col = viridis(256))
  map("world2", add = TRUE, 
      col = "grey90", lwd = border_lwd)
  
  imagePlot(x = lon, y = lat,
            rho[,,idx], main = "Rho", col = viridis(256))
  map("world2", add = TRUE, 
      col = "grey90", lwd = border_lwd)
  par(mfrow = c(1,1))
}

# STUN plots
plot_params(
  idx = 1, 
  aq_fields = df_AQ, 
  kappa2 = kappa2_s, 
  theta = theta_s, 
  rho = rho_s
)
# take log(awght - 4) to have log(range)
# -----------------------------------------------------------------------
# SETTING UP IN LATTICEKRIG
# -----------------------------------------------------------------------

chosen_day <- 21

# okay now here's how we use these parameters in LK:
# we just use the STUN ones for now 
# and let's use the first day 
kappa2_1 <- kappa2_s[,,chosen_day]
theta_1 <- theta_s[,,chosen_day]
rho_1 <- rho_s[,,chosen_day]
rhox_1 <- sqrt(rho_1)
rhoy_1 <- 1/rhox_1

# some grid setup 
rows <- dim(df_AQ)[1]
gridList<- list( x= seq( 1,rows,length.out= rows),
                 y= seq( 1,rows,length.out= rows) )
sGrid<- make.surface.grid(gridList)

# create H tensor out of params
H11 <- ( rhox_1^2 * (cos(theta_1))^2) + ( rhoy_1^2 * (sin(theta_1))^2 ) 
H12 <- (rhoy_1^2 - rhox_1^2)*(sin(theta_1)*cos(theta_1))
H21 <- H12 
H22 <- (rhox_1^2 * (sin(theta_1))^2) + (rhoy_1^2 * (cos(theta_1))^2)

# fill the high dimensional stencil (9 fields)
stencil_tensor <- array( NA, c( rows,rows,9))
stencil_tensor[,,1] <- 0.5 * H12
stencil_tensor[,,2] <- -H22
stencil_tensor[,,3] <- -0.5 * H12
stencil_tensor[,,4] <- -H11
stencil_tensor[,,5] <- kappa2_1 + 2 * H11 + 2 * H22
stencil_tensor[,,6] <- -H11
stencil_tensor[,,7] <- -0.5 * H12
stencil_tensor[,,8] <- -H22
stencil_tensor[,,9] <- 0.5 * H12

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

awght_obj <- list( x= gridList$x,  y= gridList$y, z=stencil_tensor )
class( awght_obj)<- "multivariateSurfaceGrid"


# now inside of your LKinfo, set this to be your awght object, like so:
LKinfo_test <- LKrigSetup(sGrid, NC =rows,
                     nlevel = 1, 
                     a.wghtObject =  awght_obj, 
                     normalize=FALSE, 
                     NC.buffer = 0, overlap = 2.5, nu = 1) 


# -----------------------------------------------------------------------
# PLOTTING MODEL CORRELATIONS
# -----------------------------------------------------------------------


plot_cov_surface <- function(
    LKinfo_test, 
    real_field,
    loc_x, 
    loc_y, 
    plotting = TRUE, 
    map_contour = TRUE, 
    map_lwd = 1.3, 
    cov_countour = TRUE,
    cov_lwd = 1.3
){
  cov_surface <- LKrig.cov(sGrid, rbind(c(loc_x,loc_y)), LKinfo_test)
  cov_surface <- matrix(cov_surface, nrow = rows, ncol = rows)
  
  if (plotting == TRUE){
    par(mfrow = c(1,2), mar=  c(5.1, 4.1, 4.1, 2.1))
    image.plot(x = lon, y = lat, cov_surface, col = viridis(256), 
               main = "Covariance")
    if (cov_countour == TRUE){
      contour( 
        x = lon, 
        y = lat, 
        cov_surface, 
        col = "white", 
        add = TRUE, lwd = cov_lwd
      )
    }
    if (map_contour == TRUE){
      map("world2", add = TRUE,
          col = "grey90", lwd = map_lwd)
    }
    
    imagePlot(x = lon, y = lat, real_field, 
              main = "NO2 Residual Field", col = turbo(256))
    if (cov_countour == TRUE){
      contour( 
        x = lon, 
        y = lat, 
        cov_surface, 
        col = "black", 
        add = TRUE, lwd = cov_lwd
      )
    }
    if (map_contour == TRUE){
      map("world2", add = TRUE,
          col = "grey20", lwd = map_lwd)
    }
    par(mfrow = c(1,1))
  }
  # return (cov_surface)
}


# example correlations in the po valley 
po_x <- 48
po_y <- 98
plot_cov_surface(
  LKinfo = LKinfo_test, 
  real_field = df_AQ[,,chosen_day], 
  loc_x = po_x, 
  loc_y = po_y, 
  map_contour = TRUE, 
  map_lwd = 1.3,
  cov_countour = FALSE,
  cov_lwd = 1
)

# example correlations in the ocean 

# ocean_x <- 76
# ocean_y <- 48
ocean_x <- 70
ocean_y <- 60

plot_cov_surface(
  LKinfo = LKinfo_test, 
  real_field = df_AQ[,,chosen_day], 
  loc_x = ocean_x, 
  loc_y = ocean_y, 
  map_contour = TRUE, 
  map_lwd = 1.3,
  cov_countour = TRUE,
  cov_lwd = 1
)


# random choice 
# (plot a bunch of these and see if you like the way they look)

rand_x <- sample(1:128, 1)
rand_y <- sample(1:128, 1)

plot_cov_surface(
  LKinfo = LKinfo_test, 
  real_field = df_AQ[,,chosen_day], 
  loc_x = rand_x, 
  loc_y = rand_y, 
  map_contour = TRUE, 
  map_lwd = 1.3,
  cov_countour = FALSE,
  cov_lwd = 1
)












# example of estimating an additional kappa2

# let's first take the first day residuals
NO2_field <- df_AQ[,,chosen_day]
image.plot(NO2_field, col = turbo(256))

# now let's randomly remove 80% of the field
NO2_field <- as.vector(NO2_field)
set.seed(777)
obs_idx <- sample(1:length(NO2_field), size = round(0.1 * length(NO2_field)), replace = FALSE)
NO2_obs <- NO2_field[obs_idx]
obs_locs <- sGrid[obs_idx, ]

# now let's fit normal stationary LK to it
base_model <- LatticeKrig(
  obs_locs,
  NO2_obs,
  nlevel = 1,
  NC = 128,
  normalize = FALSE,
  NC.buffer = 0,
  overlap = 2.5,
  findAwght = T
)


base_surf <- predict(base_model, sGrid)
image.plot(as.surface(sGrid, base_surf), col = turbo(256), main = "stat naive fit")

# now let's fit nonstationary using LK (no additional kappa2)
nonstat_model <- LatticeKrig(
  obs_locs,
  NO2_obs,
  LKinfo = LKinfo_test
)

nonstat_surf <- predict(nonstat_model, sGrid)
image.plot(as.surface(sGrid, nonstat_surf), col = turbo(256), main = "nonstat default fit")






# now let's do the additional kappa2 modification
LatticeKrig_hacked<- function(x, y, Z=NULL, weights=NULL,   nlevel=3, findAwght=FALSE, 
                       LKinfo=NULL, X=NULL, U=NULL, na.rm=TRUE,
                       tol=.005, verbose=FALSE, ...){
  # a crisp wrapper where many default values are exercised.
  x<- as.matrix(x)
  y<- as.matrix(y)
  if( is.null(weights)){
    weights<- rep( 1, nrow( y))
  }
  # adjust for missing values
  ind<- is.na(y)
  if( any(ind)){
    if( na.rm){
      x<- x[!ind,]
      y<- y[!ind]
      weights<- weights[!ind]
      warning("NAs removed")
      if( !is.null(Z)){
        Z<- as.matrix( Z)[!ind,]
      }
    }
    else{
      stop("NAs in y")
    }
  }
  #      	                     
  if( is.null(LKinfo) ){
    argList<-list( ...)
    # determine the geometry/dimension if not specified
    # set up some thin plate spline like default models for just Euclidean spatial domains
    # in 1,2 and 3 dimensions.              
    argList<- LatticeKrigEasyDefaults(argList,nlevel,x)
    if(verbose){
      cat("extra args:", fill=TRUE)
      print( names(argList))
    }
    LKinfo<- do.call( "LKrigSetup", c( list( x=x,  nlevel=nlevel,
                                             verbose=FALSE), argList ) )
  }  
  if( verbose){
    print(LKinfo)
  } 
  # find lambda and/ or Awght   
  if( !findAwght){
    obj<- LKrigFindLambda( x=x,y=y, weights=weights,
                           X=X, U=U, Z=Z, LKinfo=LKinfo,
                           tol=tol,
                           verbose=verbose)
    LKinfo <- LKinfoUpdate( LKinfo, lambda= obj$lambda.MLE)
  }
  else{
    obj<- LKrigFindLambdaAwght_hacked( x=x,y=y, X=X, weights=weights, U=U, Z=Z, LKinfo=LKinfo,
                                verbose=verbose)
    LKinfo <- LKinfoUpdate( LKinfo, lambda= obj$lambda.MLE,
                            a.wght=obj$a.wght.MLE)
  }   
  # final call to LKrig to get all the summary statistics 
  obj2<- c(  LKrig( x, y, weights=weights, Z=Z, X=X, U=U, LKinfo=LKinfo), 
             list(MLE= obj) )             
  class( obj2)<- c(  "LatticeKrig", "LKrig")
  obj2$call<- match.call()
  obj2$findAwght<- findAwght
  return( obj2)
}



LKrigFindLambdaAwght_hacked <- function(x, y, ...,  LKinfo,
                                 use.cholesky=NULL, 
                                 lowerBoundLogLambda =-16,
                                 upperBoundLogLambda = 4,
                                 lowerBoundOmega = -3,
                                 upperBoundOmega =  .75,
                                 factr=1e7,
                                 pgtol=1e-1,
                                 maxit=15,
                                 verbose=FALSE) {
  #require(stats)
  
  # For rectangle omega = log(kappa) = log(sqrt(Awght-4))
  # but will change with other models. 
  # Parts of the LKrig call that will be fixed.  (except updates to LKinfo)                             
  if( any( attr(LKinfo$a.wght,"isotropic") ) == FALSE  ){
    stop( paste(attr(LKinfo$a.wght,"isotropic"), 
                "findAwght only setup to estimate a single a.wght
                 parameter in the model.")
    )
  }
  LKrigArgs <- c(list(x = x, y = y), list( ...),
                 list( LKinfo = LKinfo,
                       NtrA = 0, 
                       getVarNames = FALSE)
  )
  
  if( verbose){
    cat( "LKrigFindLambdaAwght: Set of LKrigArgs before first call:", names(LKrigArgs ), fill=TRUE)
  }
  # Set up initial values of Awght and omega
  Awght.init <- as.numeric(LKrigArgs$LKinfo$a.wght[1])
  omega.start <- Awght2Omega( Awght.init, LKinfo)
  
  # Set up initial values of lambda and log lambda
  lambda.start <- LKrigArgs$LKinfo$lambda 
  
  if( is.na(lambda.start) ){ 
    llambda.start <- -1 
  } 
  else{
    llambda.start<- log( lambda.start)
  }
  # 
  if ( (llambda.start < lowerBoundLogLambda) || 
       (llambda.start > upperBoundLogLambda) || 
       (is.na(llambda.start))
  ){
    stop("Given lambda value is out of bounds.")
  }  
  #
  if(verbose){
    cat("LKrigFindLambdaAwght: llambda.start:",  llambda.start, "a.wght.start:", Awght.init, fill=TRUE)
  }
  a.wghtTemp<- omega2Awght(omega.start, LKinfo)
  lambdaTemp<- exp( llambda.start)
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
      verbose = FALSE)))
  # Update the LKrigArgs with cholesky decomp and wU  
  LKrigArgs$use.cholesky<- LKrigObject$Mc
  LKrigArgs$wU<- LKrigObject$wUb
  #
  capture.evaluations <-  rbind( c(lambdaTemp,llambda.start,
                                   a.wghtTemp, omega.start,
                                   LKrigObject$sigma2.MLE.FULL,
                                   LKrigObject$tau.MLE.FULL,
                                   LKrigObject$lnProfileLike.FULL) 
  )
  if(verbose){
    cat("Capture.evaluations first call", fill=TRUE )
    cat("lambda", "log lambda", "a.wght", "omega",
        "sigma2MLE", "tauMLE", "logProfileLike", fill=TRUE)
    cat( capture.evaluations, fill=TRUE)
  }
  
  #####  optimze likelihood over log lambda  and over omega =  log( a.wght -4)/2
  capture.env <- environment()
  # last two arguments are specific to this objectinve function     
  result <- try(optim(c(llambda.start, omega.start),
                      LambdaAwghtObjectiveFunction_hacked, 
                      lower=c(lowerBoundLogLambda,lowerBoundOmega), 
                      upper=c(upperBoundLogLambda,upperBoundOmega), 
                      method="L-BFGS-B",
                      #                      method="BFGS",
                      control=list(fnscale = -1,factr=factr,
                                   pgtol=pgtol, maxit=maxit,
                                   ndeps = c(.05,.05)),
                      
                      LKrigArgs=LKrigArgs,
                      capture.env= capture.env,
                      verbose=verbose
  ))
  if(verbose){
    cat("Results from optimize:", fill=TRUE)
    print( result )
  }
  evalSummary <- !(class( result)== "try-error")
  llambda.MLE <- result$par[1]
  lambda.MLE<- exp( llambda.MLE)
  omega.MLE<- result$par[2]
  a.wght.MLE<- omega2Awght(omega.MLE,  LKrigArgs$LKinfo )
  LKrigArgs$NtrA <- 20
  
  LKrigObject <- do.call("LKrig_hacked", c(LKrigArgs,
                                    list(
                                      lambda = lambda.MLE,
                                      a.wght = a.wght.MLE
                                    )
  )
  )
  
  ###### end optimze block    
  # save summary results from this set of parameters.
  # Output to be saved     
  out <- rep(NA, 11)
  names( out) <-  c("EffDf", "lnProfLike", "GCV", "tau.MLE", "sigma2.MLE", 
                    "lambda.MLE", "a.wght.MLE", "lnLike", "functionEval", 
                    "gradientEval", "totalEval")
  out[ 1] <- LKrigObject$trA.est
  out[ 2] <- LKrigObject$lnProfileLike.FULL
  out[ 3] <- LKrigObject$GCV
  out[ 4] <- LKrigObject$tau.MLE.FULL
  out[ 5] <- LKrigObject$sigma2.MLE.FULL
  out[ 6] <- lambda.MLE
  out[ 7] <- a.wght.MLE
  out[ 8] <- LKrigObject$lnLike.FULL
  out[ 9] <- result$counts[1]
  out[10] <- result$counts[2]
  out[11] <- nrow(capture.evaluations )
  
  # Name columns  of likelihood eval. 
  dimnames(capture.evaluations)<- list( NULL,
                                        c("lambda","logLambda","a.wght","omega",
                                          "sigma2.MLE", "tau.MLE","lnProfileLike.FULL"))
  return(list(summary = out,
              LKinfo = LKrigObject$LKinfo,
              llambda.start = llambda.start,
              Awght.start=omega2Awght( omega.start,LKrigArgs$LKinfo),
              lambda.MLE = lambda.MLE,
              a.wght.MLE = a.wght.MLE,
              omega.MLE =omega.MLE,
              llambda.MLE=llambda.MLE,
              lnLike.eval = capture.evaluations,
              call = match.call() )
  )        
}

# Define the objective function 
LambdaAwghtObjectiveFunction_hacked<- function(PARS, LKrigArgs, capture.env, verbose=FALSE ) {
  lambdaTemp <- exp( PARS[1] )
  a.wghtTemp <-  omega2Awght( PARS[2], LKrigArgs$LKinfo)
  
  hold <- do.call("LKrig_hacked",          
                  c(LKrigArgs, list(
                    lambda = lambdaTemp,
                    a.wght = a.wghtTemp)
                  )
  )[c("sigma2.MLE.FULL","tau.MLE.FULL","lnProfileLike.FULL")] 
  rowForCapture <-c( lambdaTemp,PARS[1],
                     a.wghtTemp,PARS[2],
                     hold$sigma2.MLE.FULL,
                     hold$tau.MLE.FULL,
                     hold$lnProfileLike.FULL 
  )
  if( verbose){
    cat( rowForCapture, fill=TRUE )
  }
  lnProfileLike.FULL<- hold$lnProfileLike.FULL 
  temp.eval <- get("capture.evaluations",
                   envir = capture.env )
  assign("capture.evaluations",rbind(temp.eval, rowForCapture),
         envir = capture.env)
  return(lnProfileLike.FULL)
}



LKrig_hacked <- function(x, y,
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
  if ( is.null(LKinfo) ) {
    LKinfo <- do.call("LKrigSetup", c(list(x = x ), list(...),
                                      list(verbose = verbose)))
  }
  else{
    LKinfo<- do.call("LKinfoUpdate", c(list(LKinfo=LKinfo), list(...)) )	
  }
  
  if( verbose){
    cat(" ", fill=TRUE)
    cat("LKrig: updated LKinfo object", fill=TRUE)
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
  if( getVarNames){
    xName<- as.character( substitute( x) )
    ZName<- as.character( substitute( Z) )
    UName<- as.character( substitute( U) )
    xName<- tail( xName, 1)
    ZName<- tail( ZName, 1)
    UName<- tail( UName, 1)
    
    # just take last component
    
  }
  else{
    xName<- "xVar"
    ZName<- "ZVar"
    UName<- "UVar"
  }  
  #
  # create the initial parts of LKrig object
  # this list is added to as the computation proceeds 
  # using the device  object<- list( object, newStuff)
  # and the full object is only obtained at the end 
  # NOTE default for weights are just 1's and filled in by 
  # the next call    
  
  object<- createLKrigObject( x, y, 
                              weights = weights,
                              Z = Z,
                              X = X,
                              U = U, 
                              LKinfo = LKinfo,
                              xName = xName, 
                              ZName = ZName,
                              UName = UName, 
                              verbose = verbose)
  
  nObs <-  nrow( object$y )
  nReps <- ncol( object$y )
  # for readablity make a local copy of LKinfo
  # but don't change it in this function! 
  LKinfo<- object$LKinfo                                              	 	
  # Begin computations ....
  # weighted observation vector
  wy <- sqrt(object$weights) * object$y
  
  # create matrix for fixed part of model    
  # Spatial drift matrix -- default is assumed to be linear in coordinates (m=2)
  # and includes possible covariate(s) -- the Z matrix.
  # the choice of fixed part of the model is controlled in LKinfo
  # (see also LKrigSetup)
  if (is.null(wU)) {
    wU<- LKrigMakewU( object,  verbose=verbose)
  }
  
  # some column indices to keep track of fixed part of the model	
  # NOTE nZ <= nt because Z is a subset of U
  object$nt <- ifelse( is.null(ncol(wU)), 0, ncol(wU))
  # create matrix for random part of model (basis functions)
  #  wX is the matrix of sum( N1*N2) basis function (columns) evaluated at the N locations (rows)
  # and multiplied by square root of diagonal weight matrix
  # this can be a large matrix if not encoded in sparse format.
  if (is.null(wX)) {
    timewX<- system.time(
      wX<- LKrigMakewX( object, verbose=verbose)
    )	
  }
  else{
    timewX<- rep(0,5)
  }	
  
  #   Precision matrix of the lattice process
  #   inverse of Q is proportional to the covariance matrix of the Markov Random Field
  timeQ<-system.time(
    # THIS IS HONESTLY THE ONLY BIG CHANGE I THINK -----------------------------------------
    # Q <- LKrig.precision(LKinfo, verbose=verbose)
    Q <- nonstat_Q
  )
  print(LKinfo$a.wght[[1]]-4)
  # print(diag(Q))
  diag(Q) <- diag(Q) + (LKinfo$a.wght[[1]]-4)
  # ENDS HERE ------------------------------------------------------------------------------
  
  if( LKinfo$dense){
    if( !is.null(  use.cholesky)){
      stop("Can not update (use.cholesky) with dense matrix option")
    }
    Q<- spam2full(Q)
    wX<- spam2full(wX)
  }
  
  if( verbose & (!LKinfo$dense) ) {
    cat("LKrig: Nonzero entries in Q:", length(Q@entries), fill=TRUE)		
  }
  
  # G is the regularized (ridge) regression matrix that is 
  # the key to the entire algorithm:
  timeM<- system.time(	
    G <- t(wX) %*% wX + LKinfo$lambda * (Q)
  )
  if( verbose){
    if( !LKinfo$dense){
      cat("LKrig: Nonzero entries in M:", length(G@entries), fill=TRUE)	
    }
    else{
      cat( "Dense matrix methods used", fill=TRUE)
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
    timeChol<- system.time(
      GCholesky <- chol(G, memory = LKinfo$choleskyMemory)
    )
  } else {	
    timeChol<- system.time(
      GCholesky <- update.spam.chol.NgPeyton(use.cholesky, G)
    )
  }
  
  if( !LKinfo$dense){
    nonzero.entries<- length(GCholesky@entries)
  }
  else{
    nonzero.entries<-NA
  }
  
  
  if( verbose){
    cat("LKrig: nonzero entries of GCholesky:",nonzero.entries, fill=TRUE)
  }
  # use GCholesky to find coefficients of estimate
  # Note that this functions also finds an important piece of the likelihood (quad.form)
  timeCoef<- system.time(
    out1 <- LKrig.coef(GCholesky, wX, wU, wy,
                       LKinfo$lambda, 
                       collapseFixedEffect = LKinfo$collapseFixedEffect, 
                       verbose=verbose)
  )
  
  # Note collapseFixedEffect added as component here in the return
  # finding coefficients
  # fill in names of the fixed coefficients 
  # fill in names of fixed model coefficients
  
  rownames(out1$d.coef)<- colnames( wU )
  #
  object <- c(object, out1)
  if( verbose){
    cat("fixed model coefficients", fill=TRUE)
    cat( object$d.coef, fill=TRUE)
  }
  
  # compute predicted values  and residuals
  wfitted.values <- (wX %*% out1$c.coef)
  if ( !is.null(wU) ) {
    wfitted.values.fixed <- (wU %*% out1$d.coef)
    wfitted.values <- wfitted.values.fixed + wfitted.values
  }
  # X and U actully include the weights so need to divide these
  # out to get fitted values	
  object$fitted.values<- wfitted.values/sqrt(object$weights)		
  # For reference: fitted.values <- predict.LKrig(object, x, Znew = object$Z)
  # but at this point it is less efficient because X will be recomputed.
  object$residuals <- object$y - object$fitted.values	
  # find likelihood
  timeLike<- system.time(	
    #	out2 <- LKrig.lnPlikeOLD(GCholesky, Q, wy,
    #	             object$residuals, object$weights,
    #	             LKinfo)
    out2 <- LKrig.lnPlike(GCholesky, Q, object$quad.form,
                          nObs, nReps,
                          object$weights,LKinfo)
  )
  if( verbose){
    cat("Likelihood/MLE list:",  fill=TRUE)
    print( out2)
  }   
  object <- c(object, out2)
  
  # estimate trace of hat matrix (effective degrees of freedom)
  # by Monte Carlo if NtrA greater than zero
  timeTrA<- system.time(
    if (NtrA > 0) {
      out3 <- LKrig.traceA(GCholesky, wX, wU, LKinfo$lambda, object$weights, NtrA, iseed = iseed)
      # find GCV using this trace estimate
      n<- length( object$weights)
      out3$GCV = (sum(object$weights * (object$residuals)^2)/n)/(1 - out3$trA.est/n)^2
    } else {
      out3 <- list(trA.est = NA, trA.SE = NA, GCV = NA)
    }
  )
  object <- c(object, out3)
  
  # create the table of times for individual function calls
  timingTable<- rbind(timewX, timeQ, timeM, timeChol, timeCoef, timeLike,  timeTrA)
  timingTable<- timingTable[,1:3]
  timingTable <- rbind( timingTable, colSums(timingTable))
  
  # last of required arguments to LKrig object
  object <- c(object,
              list( 
                lambda.fixed = LKinfo$lambda, 
                nonzero.entries = nonzero.entries,
                call = match.call(), 
                timingLKrig = timingTable )
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

nonstat_Q <- LKrig.precision(LKinfo_test)

nonstat_model_hacked <- LatticeKrig_hacked(
  obs_locs,
  NO2_obs,
  nlevel = 1,
  NC = 128,
  normalize = FALSE,
  NC.buffer = 0,
  overlap = 2.5,
  findAwght = T
)

nonstat_surf_hacked <- predict(nonstat_model_hacked, sGrid)
# image.plot(as.surface(sGrid, nonstat_surf_hacked), col = turbo(256))

print(nonstat_model_hacked$LKinfo$a.wght)



length(LKinfo_test$a.wght[[1]][,5])

LKinfo_test$a.wght[[1]][,5] <- LKinfo_test$a.wght[[1]][,5] + (nonstat_model_hacked$LKinfo$a.wght[[1]]-4)

nonstat_model_hacked_add  <- LatticeKrig(
  obs_locs,
  NO2_obs,
  LKinfo = LKinfo_test
)

nonstat_surf_hacked_add <- predict(nonstat_model_hacked_add, sGrid)
image.plot(as.surface(sGrid, nonstat_surf_hacked_add), col = turbo(256),
           main = "nonstat hacked fit")

image.plot(as.surface(sGrid, nonstat_surf_hacked_add - nonstat_surf), col = turbo(256),
           main = "diff")

nonstat_model$LKinfo$a.wght[[1]][,5] - nonstat_model_hacked_add$LKinfo$a.wght[[1]][,5]

nonstat_model_hacked$LKinfo$a.wght[[1]]
base_model$LKinfo$a.wght[[1]]
mean(4 + kappa2_s[,,chosen_day])
median(4 + kappa2_s[,,chosen_day])

