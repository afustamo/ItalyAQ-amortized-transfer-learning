
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
  if (verbose){
    print(LKinfo$a.wght[[1]]-4)
  }
  # print(diag(Q))
  diag(Q) <- diag(Q) + (kappa2_weights * (LKinfo$a.wght[[1]]-4))
  # print(kappa2_weights)
  # print(LKinfo$a.wght[[1]]-4) # for debugging
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
  # timeM<- system.time(	
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