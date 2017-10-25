# Implements the ChiSq backend. Currently all the hevy lifting is done via the bmdupdate library.
#library(bmdupdate)

#' @importFrom methods new 
backend.chisq <- function(X, Y, parallel = FALSE, ...) {
  p <- backend.base(X, Y, ...)
  p$obj = new(bmdupdate::BmdUpdater, X, Y)
  p$two_sided = TRUE
  p$parallel <- parallel
  class(p) <- c("chisq", class(p))
  return(p)
}

pvals.chisq <- function(p, B){
  p$obj$pvals(B, p$parallel)
}

# pvalsR_chisq <- function (B) {
#   
#   if (length(B) == 0)
#     return(integer(0))
#   
#   test_X <- min(B) > dx
#   nFixd <- length(B)
#   
#   if (test_X) {
#     
#     # Getting fixed matrix
#     fixdIndx <- match(B, Yindx)
#     fixdMat <- Y[ , fixdIndx, drop = FALSE]
#     
#     if (nFixd > n) {
#       
#       # Do n^2|B| computation
#       
#     } else {
#       
#       # Do n|B|^2 computation
#       
#     }
#     
#   } else {
#     
#     # Getting indices
#     fixdIndx <- match(B, Xindx)
#     fixdMat <- X[ , fixdIndx, drop = FALSE]
#     
#     if (nFixd > n) {
#       
#       # Do n^2|B| computation
#       
#     } else {
#       
#       # Do n|B|^2 computation
#       
#     }
#     
#   }
#   
#   # Compute as and bs
#   corsums <- as.vector(rowSums(xyCors))
#   zstats <- sqrt(n) * corsums / sqrt(allvars)
#   if (twoSided) {
#     pvals <- 2 * pnorm(abs(zstats), lower.tail = FALSE)
#   } else {
#     pvals <- pnorm(zstats, lower.tail = FALSE)
#   }
#   return(pvals)
#   
# }