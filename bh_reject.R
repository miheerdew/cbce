bh_rejectR <- function (pvals, alpha, conserv = TRUE) {
  
  m <- length(pvals)
  
  if (!conserv) {
    pvals_adj <- m * pvals / rank(pvals, ties.method = "first")
  } else {
    mults <- sum(1 / c(1:m))
    pvals_adj <- mults * m * pvals / rank(pvals, ties.method = "first")
  }
  
  if (sum(pvals_adj <= alpha) > 0) {
    thres <- max(pvals[pvals_adj <= alpha])
    return(which(pvals <= thres))
  } else {
    return(integer(0))
  }
  
}