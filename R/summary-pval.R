summary_perm_moments <- function(X, Y) {
  
  n <- nrow(X)
  
  X <- scale(X)/sqrt(n-1)
  Y <- scale(Y)/sqrt(n-1)
  
  lambda1 <- symmetric_eigenvals(X)
  lambda2 <- symmetric_eigenvals(Y)
  
  # Equal to sum(lambda) since X, Y were normalized.
  sum.1 <- ncol(X)
  sum.2 <- ncol(Y)
  
  sum.sq.1 <- sum(lambda1^2)
  sum.sq.2 <- sum(lambda2^2)
  
  sum.cube.1 <- sum(lambda1^3)
  sum.cube.2 <- sum(lambda2^3)
  
  #Degree of freedom
  df <- nrow(X) - 1
  
  #Central moments.
  mu1 <- sum.1 * sum.2 / df
  
  a <- 3 - (df-3)/(df-1)
  b <- (df+1)/(df-1)
  c <- 2/(df-1)
  
  mu2 <- (a*sum.sq.1*sum.sq.2 +
            b*(sum.1*sum.2)^2 -
            c*(sum.sq.1*sum.2^2 + sum.sq.2*sum.1^2)
  )/(df*(df+2))
  
  
  # --------- Preparation for the third moment -------
  a1 <- 1/(df*(df+2)*(df+4))
  a2 <- (df+3)/((df-1)*df*(df+2)*(df+4))
  a3 <- ((df^2 + 3*df -2)/((df+4)*(df+2)))/(df*(df-1)*(df-2))
  
  a <- a3*sum.1^3 +
    3*(a2 - a3)*(sum.sq.1)*(sum.1) +
    (a1 - a3 - 3*(a2 - a3))*sum.cube.1
  
  b1 <- 3/(df*(df+2)*(df+4))
  b2 <- 3*(df+3)/((df-1)*df*(df+2)*(df+4))
  b3 <- (df+1)/((df-1)*df*(df+2)*(df+4))
  b4 <- (df+3)/((df-1)*df*(df+2)*(df+4))
  
  b <- (b1 - b2 - 2*b3 + 2*b4)*sum.cube.1 +
    (b2 + 2*b3 - 3*b4)*sum.sq.1*sum.1 +
    b4*sum.1^3
  
  c1 <- 15/(df*(df+2)*(df+4))
  c2 <- 3/(df*(df+2)*(df+4))
  c3 <- 1/(df*(df+2)*(df+4))
  
  c <- c3*sum.1^3 +
    (c1 - c3 - 3*(c2 - c3))*sum.cube.1 +
    3*(c2-c3)*sum.sq.1*sum.1
  
  # Finally, the third moment:
  mu3 <- a*sum.2^3 +
    (c - a - 3*(b-a))*sum.cube.2 +
    3*(b-a)*sum.sq.2*sum.2
  
  mu <-  mu1 #Mean
  var <- mu2 - mu1^2
  skewness <- (mu3-3*mu*var-mu^3)/(var^1.5) #skewness
  
  c(mean=mu, var=var, skewness=skewness)
}

symmetric_eigenvals <- function(X) {
  m <- ncol(X)
  n <- nrow(X)
  A <- if(m <= n) crossprod(X) else tcrossprod(X)
  eigen(A)$values
}

#' Compute the score (i.e a summarizing p-vlaue) for the 
#' bimodule given by the entire X and Y matrices.
#' 
#' @export
#' @keywords internal
summary_pval <- function(X, Y, log.p = TRUE, summary=NULL) {
  pm <- as.list(summary_perm_moments(X, Y))
  
  b <- 8/pm$skewness^2
  a <- sqrt(pm$var/(2*b))
  d <- pm$mean - a*b
  
  if(is.null(summary)) {
    summary <- sum(cor(X, Y)^2)
  }
  
  pchisq((summary-d)/a, df=b, lower.tail=FALSE, log.p = log.p)
}
