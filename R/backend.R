# backends should implement pvals, cors and pvals_singleton.
# see the base backend which implements cors and pvals_singleton and can be used by other backends

#' Calculate p-values against a testing set.
#' 
#' Given a testing set B, return the p-values for the opposite side.
#' This is an S3 method that is implemented by all backends.
#' 
#' @param bk An object of class backend
#' @param B A set of either X or Y indices (using global numbering). The method calculates the pvalues from this set to the opposite side. 
#' 
#' @return The result is a vector of p-values of length ncol(X) or ncol(Y) depending on whether B is from Y or X respectively.
#' @seealso \code{\link{pvals_singleton}}
#'@export
pvals <- function(bk, B) {
  if (length(B) == 0) return(integer(0))
  UseMethod("pvals", bk)
}

#' Calculate p-values against a single variable.
#' 
#' Given a node return the p-values for the opposite side for the correlations to this node. This is supposed to be a special case of \code{\link{pvals}} (with a singleton set) but it is more faster and accurate.
#' 
#' @param bk An object of class backend
#' @param indx Either an X or Y indices (using global numbering).  
#' @return The result is a vector of p-values of length ncol(X) or ncol(Y) depending on whether indx is from Y or X respectively.
#' @seealso \code{\link{pvals}}
#'@export
pvals_singleton <- function(bk, indx) {
  #Given indx return the set of all p-values to opposite side. 
  #This is certainly a special case of pvals(b, A) with A = indx, but can be computed much faster since it does not involve sums of correlations.
  UseMethod("pvals_singleton", bk)
}

cors <- function(bk, A) {
  #calculate the correlations from set A.``
  UseMethod("cors", bk)
}

init <- function(p, indx, alpha, init_method) {
  pvals <- pvals_singleton(p, indx) 
  switch(init_method,
        "conservative-BH" = bh_reject(pvals, alpha, conserv = TRUE),
        "non-conservative-BH" = bh_reject(pvals, alpha, conserv = FALSE),
        "BH-0.5" =  bh_reject(pvals, 0.5, conserv = TRUE),
        "no-multiple-testing" = which(pvals <= alpha),
        stop(paste("Unknown init_method:", init_method))
  )
}

## Defaults
pvals_singleton.default <- function(bk, indx) {
  pvals(bk, indx)
}