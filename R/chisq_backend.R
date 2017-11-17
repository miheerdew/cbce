# Implements the ChiSq backend. Currently all the hevy lifting is done via the bmdupdate library.

#' @describeIn backend Constructor for sum of squared correlations under uncorrelated Gene ans SNP sets (the weak null).
#' @param parallel Logical Run pvalue computation using multiple threads?
#' @inheritParams backend.base
#'@importFrom methods new 
#'@export
backend.chisq <- function(X, Y, parallel = FALSE, calc_full_cor=FALSE) {
  p <- backend.base(X, Y, calc_full_cor)
  p$obj = new(bmdupdate::BmdUpdater, X, Y)
  p$two_sided = TRUE
  p$parallel <- parallel
  class(p) <- c("chisq", class(p))
  return(p)
}

#' @describeIn pvals C implementation for sum of squared correlations under uncorrelated Gene and SNP sets (the weak null)
#' @export
pvals.chisq <- function(bk, B){
  bk$obj$pvals(B, bk$parallel)
}