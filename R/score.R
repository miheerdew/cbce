score <- function(X, Y) {
  X <- scale(X)
  Y <- scale(Y)
  dx <- ncol(X)
  dy <- ncol(Y)
  n <- nrow(X)
  N <- n - 1
  
  X2 <- X^2
  Y2 <- Y^2
  
  XYs <- do.call("cbind", rlist::list.map(1:dy, X*Y[,.]))
  X2s <- matrix(rep(X2, times=dy), nrow=n)
  Y2s <- do.call("cbind", rep(rlist::list.map(1:dy, Y2[,.]), each=dx))
  
  #The indexing is (i, j) -> i + (j-1)*dx
  Rs <- as.vector(cor(X, Y))
  
  S1 <- crossprod(XYs)/N
  S2 <- tcrossprod(Rs) * crossprod(X2s + Y2s)/N
  S3 <- crossprod(X2s + Y2s, XYs)/N
  
  S <- S1 + S2/4 - (S3 + t(S3))/2
  CompQuadForm::davies(n*sum(Rs^2), base::eigen(S)$values, sigma=1, acc=1e-7)
}
