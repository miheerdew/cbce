# Implements the Normal backend in C. All the hevy lifting is done via the bmdupdate library.

#' @describeIn backend Constructor for sum of correlations under uncorrelated Gene ans SNP sets (the weak null).
#' @param parallel Logical Run pvalue computation using multiple threads?
#' @inheritParams backend.base
#'@importFrom methods new
#'@export
backend.normalc <- function(X, Y, parallel = FALSE, calc_full_cor=FALSE) {
  p <- backend.base(X, Y, calc_full_cor)
  p$obj = new(bmdupdate::BackendNormal, X, Y)
  p$two_sided = FALSE
  p$parallel <- parallel
  p$normal_vector_pval = function(v) stats::pnorm(sum(v), sd=sqrt(length(v)), lower.tail=FALSE)
  class(p) <- c("normalc", class(p))
  return(p)
}

#' @describeIn pvals C implementation for sum of correlations under uncorrelated Gene and SNP sets (the weak null)
#' @export
pvals.normalc <- function(bk, B) {
    bk$obj$pvals(B, bk$parallel)
}

#' @describeIn pvals C implementation for sum of correlations under uncorrelated Gene and SNP sets (the weak null)
#' @export
pvals_quick.normalc <- function(bk, B) {
  R <- cors(bk, B)
  stats::pnorm(sqrt(bk$n) * rowSums(R), sd=sqrt(length(B)), lower.tail = FALSE)
}
