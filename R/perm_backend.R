# Implements the Perm backend using Fred's moment based approximation.

#' @describeIn backend Constructor for sum of squared correlations under independent Gene and SNP sets (the strong null).
#' @inheritParams backend.base
#'@importFrom methods new
#'@export
backend.perm <- function(X, Y, calc_full_cor=FALSE) {
  p <- backend.base(X, Y, calc_full_cor)
  p$two_sided <- TRUE
  # No need to change:
  #  p$normal_vector_pval
  class(p) <- c("perm", class(p))
  return(p)
}

#' @describeIn pvals Implementation for sum of squared correlations under independent Gene and SNP sets. 
#'             Uses the permutation moments approximation.
#' @export
pvals.perm <- function(bk, B, thresh.alpha=1) {
  m <- length(B)
  n <- bk$n

  X <- (if(min(B) > bk$dx) bk$Y[,B-bk$dx] else bk$X[,B])/sqrt(n-1)
  A <- if(m <= n) crossprod(X) else tcrossprod(X)
  lambda <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  
  # The code for the moments of the permutation distribution.
  mu <- m/(n-1) #Mean
  sigma2 <- (sum(lambda^2) - m^2/(n-1))*(2/(n^2-1)) #Variance
  c <- 1/((n^2-1)*(n+3))
  mu3 <- c*(m^3+6*m*sum(lambda^2)+8*sum(lambda^3)) #non-central 3rd moment
  gamma1 <- (mu3-3*mu*sigma2-mu^3)/(sigma2^1.5) #skewness
  
  # The code to fit a shifted Chi-Square
  # r^2 ~ a + d \chi_b^2
  b<-8/gamma1^2
  a<-sqrt(sigma2/(2*b))
  d<-mu-a*b
  
  R <- cors(bk, B)
  Tstat <- (rowSums(R^2) - d)/a
  
  pvals <- rep(NA, length(Tstat))
  thresh <- qchisq(thresh.alpha, df=b, lower.tail=FALSE)
  
  interesting <- Tstat > thresh
  
  pvals[interesting] <- pchisq(Tstat[interesting], df=b, lower.tail=FALSE)
  
  return(pvals)
}

pvals_quick.perm <- function(bk, B) {
  R <- cors(bk, B)
  stats::pchisq(bk$n * rowSums(R^2), length(B), lower.tail = FALSE)
  
}

#' @describeIn pvals_singleton implementation for the perm class
#' @export
pvals_singleton.perm <- function(bk, indx, thresh.alpha=1) {
  # An easy way to calculate p-values from an indx
  a <- 0.5
  b <- 0.5*(bk$n - 2)
  
  thresh <- qbeta(thresh.alpha, a, b, lower.tail = FALSE)

  r2 <- cors(bk, indx)^2
  interesting <- r2 > thresh

  pvals <- rep(NA, length(r2))
  pvals[interesting] <- pbeta(r2[interesting], a, b, lower.tail = FALSE)
  
  return(pvals)
}


