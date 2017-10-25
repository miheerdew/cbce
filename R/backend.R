# backends should implement pvals, cors and pvals_singleton.
# see the base backend which implements cors and pvals_singleton and can be used by other backends

pvals <- function(bk, B) {
  #Given a testing set B, return the p-values for the opposite side
  if (length(B) == 0) return(integer(0))
  UseMethod("pvals", bk)
}

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