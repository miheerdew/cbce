# Implements the Naive (under Gene/SNP independence) versions of p-value computations

#'@describeIn backend Constructor for sum of correlations (one sided) under the strong independece (within Gene and Snips and between them).
#'@inheritParams backend.base
#'@export
backend.indepNormal <- function(X, Y, calc_full_cor=FALSE) {
  # TODO: change precomp$normal_vector_pval
  p <- backend.base(X, Y, calc_full_cor)
  p$two_sided = FALSE
  class(p) <- c("indepNormal", class(p))
  return(p)
}

#' @describeIn pvals p-values for the sum of correlations under strong independence (with in Gene and SNP's and between them too).
#' @export
pvals.indepNormal <- function(bk, B){
  R <- cors(bk, B)
  stats::pnorm(sqrt(bk$n) * rowSums(R), sd=sqrt(length(B)), lower.tail = FALSE)
}

#'@describeIn backend Constructor for sum of squared correlations under strong independence.
#'@inheritParams backend.base
#'@export
backend.indepChiSq <- function(X, Y, calc_full_cor=FALSE) {
  # Don't change p$normal_vector_pval
  p <- backend.base(X, Y, calc_full_cor)
  p$two_sided = TRUE
  class(p) <- c("indepChiSq", class(p))
  return(p)
}

#' @describeIn pvals p-values for the sum of squared correlations under strong independence (with in Gene and SNP's and between them too).
#' @export
pvals.indepChiSq <- function(bk, B){
  R <- cors(bk, B)
  stats::pchisq(bk$n * rowSums(R^2), length(B), lower.tail = FALSE)
}
