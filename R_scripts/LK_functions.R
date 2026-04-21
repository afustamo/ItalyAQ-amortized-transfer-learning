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
est_additional_kappa2 <- function(s, y, Z, LKinfo, kappa2_weights, verbose = FALSE) {
  kappa2_weights <- kappa2_weights
  
  LKinfo_additional <- LKinfo
  LKinfo_additional$a.wght <- NA
  LKinfo_additional$a.wghtObject <- NULL
  LKinfo_additional$a.wght <- list()
  LKinfo_additional$a.wght[[1]] <- 4.01
  attr(LKinfo_additional$a.wght, "isotropic") <- TRUE
  
  Lmodel_additional <- LatticeKrig_hacked(s, y, Z, LKinfo = LKinfo_additional, findAwght = T)
  
  LKinfo_new <- LKinfo
  LKinfo_new$a.wght[[1]][, 5] <- LKinfo_new$a.wght[[1]][, 5] +
    (Lmodel_additional$LKinfo$a.wght[[1]] - 4)
  
  return(LKinfo_new)
}
