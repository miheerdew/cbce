# Implements the ChiSq backend using Fred's moment based approximation.

#' @describeIn backend Constructor for sum of squared correlations under independent Gene and SNP sets (the strong null).
#' @inheritParams backend.base
#'@importFrom methods new
#'@export
backend.chisq <- function(X, Y, calc_full_cor=FALSE) {
  p <- backend.base(X, Y, calc_full_cor)
  # No need to change:
  #  p$normal_vector_pval
  class(p) <- c("chisq", class(p))
  return(p)
}

#' @describeIn pvals Implementation for sum of squared correlations under independent Gene and SNP sets. 
#'             Uses the permutation moments approximation.
#' @export
pvals.chisq <- function(bk, B) {
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
  b<-8/gamma1^2
  a<-sqrt(sigma2/(2*b))
  d<-mu-a*b
  
  R <- apply(cors(bk, B), 1, crossprod)
  pchisq((R-d)/a, df=b, lower.tail=FALSE)
}

pvals_quick.chisq <- function(bk, B) {
  R <- cors(bk, B)
  stats::pchisq(bk$n * rowSums(R^2), length(B), lower.tail = FALSE)
  
}

