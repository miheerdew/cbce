# Implements the ChiSq backend. Currently all the hevy lifting is done via the bmdupdate library.

#' @describeIn backend Constructor for sum of squared correlations under uncorrelated Gene ans SNP sets (the weak null).
#' @param parallel Logical Run pvalue computation using multiple threads?
#' @param fast_approx Use the BigX_fast approximation?
#' @inheritParams backend.base
#'@importFrom methods new
#'@export
backend.chisq <- function(X, Y, parallel = FALSE, calc_full_cor=FALSE, fast_approx=FALSE) {
  p <- backend.base(X, Y, calc_full_cor)
  p$obj = new(bmdupdate::BmdUpdater, X, Y)
  p$two_sided = TRUE
  p$parallel <- parallel
  p$fast_approx <- fast_approx
  # No need to change:
  #  p$normal_vector_pval
  class(p) <- c("chisq", class(p))
  return(p)
}

#' @describeIn pvals C implementation for sum of squared correlations under uncorrelated Gene and SNP sets (the weak null)
#' @export
pvals.chisq <- function(bk, B) {
  if(bk$fast_approx){
    bk$obj$pvals(B, bk$parallel, bmdupdate::Method$BigX_fast)
  } else {
    bk$obj$pvals(B, bk$parallel, bmdupdate::Method$Auto)
  }
}

pvals_quick.chisq <- function(bk, B) {
  R <- cors(bk, B)
  stats::pchisq(bk$n * rowSums(R^2), length(B), lower.tail = FALSE)
}

#' @describeIn mask C implementation for sum of squared correlations under uncorrelated Gene and SNP sets (the weak null)
#' @export
mask.chisq <- function(bk, Bx, By) {
  bk$obj$mask(c(Bx,By))
  mask.base(bk, Bx, By)
}
