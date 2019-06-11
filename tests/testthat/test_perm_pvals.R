context("Test perm pvals")

freds_pval <- function(B0, e=parent.frame()) {
  test <- tstat_pairs(B0, e)
  mx <- perm_moments(e$X[, B0$x])
  my <- perm_moments(e$Y[, B0$y])
  list(x=fit_shifted_chisq(my, x=test$x), y=fit_shifted_chisq(mx, x=test$y))
}

tstat_pairs <- function(B0, e=parent.frame()) {
  tx <- rowSums(cor(e$X, e$Y[, B0$y])^2)
  ty <- rowSums(cor(e$Y, e$X[, B0$x])^2)
  list(x=tx, y=ty)
}

perm_moments <- function(X) {
  #Fred's code for the moments of the permutation distribution
  m <- ncol(X)
  n <- nrow(X)
  X <- scale(X)/sqrt(n-1)
  A <- if(m <= n) crossprod(X) else tcrossprod(X)
  lambda <- eigen(A)$values
  
  mu <- m/(n-1) #Mean
  sigma2 <- (sum(lambda^2) - m^2/(n-1))*(2/(n^2-1)) #Variance
  c <- 1/((n^2-1)*(n+3))
  mu3 <- c*(m^3+6*m*sum(lambda^2)+8*sum(lambda^3)) #non-central 3rd moment
  gamma1 <- (mu3-3*mu*sigma2-mu^3)/(sigma2^1.5) #skewness
  
  c(mean=mu, var=sigma2, skewness=gamma1, mu3=mu3)
}

fit_shifted_chisq <- function(perm_moments, type='p', x=NA, nsim=NA) {
  pm <- as.list(perm_moments)
  
  # a shifted chisquare to obtain a distributional approximation
  b<-8/pm$skewness^2
  a<-sqrt(pm$var/(2*b))
  d<-pm$mean-a*b
  
  switch(type,
         d=dchisq((x-d)/a, df=b)/a,
         p=pchisq((x-d)/a, df=b, lower.tail=FALSE),
         r=rchisq(nsim, df=b)*a + d)
}

set.seed(12345)

test_that("pvals are same as simulation", {
  set.seed(12345)
  
  dx <- 5
  dy <- 6
  n <- 4
  
  X <- matrix(rnorm(dx*n), nrow=n, ncol=dx)
  Y <- matrix(rnorm(dy*n), nrow=n, ncol=dy)

  #Induce some correlation.  
  X[, 1:2] <- rowSums(X)
  Y[, 1] <- rowSums(Y)
  
  bk <- backend.perm(X, Y)
  py <- pvals(bk, 1:dx)
  px <- pvals(bk, (1:dy) + dx)
  
  e <- new.env()
  e$X <- X
  e$Y <- Y
  
  fp <- freds_pval(list(x=1:dx, y=1:dy), e)
  expect_equal(fp$x, px)
  expect_equal(fp$y, py)
})

test_df_correction <- function(dx, dy, n, ncov) {
  alpha <- 0.01
  
  z <- matrix(rnorm(ncov*n), nrow=n, ncol=ncov)
  
  X <- matrix(rnorm(dx*n), nrow=n, ncol=dx)
  Y <- matrix(rnorm(dy*n), nrow=n, ncol=dy)
  
  X[, 1:2] <- X[, 1:2] + rowSums(X)
  Y[, 1] <- Y[,1] + rowSums(Y)
  
  #Residualize the X and Y matrices w.r.t Z
  z <- scale(z, scale=FALSE)
  q <- qr(z)
  Q <- qr.Q(q)
  X <- X - Q %*% crossprod(Q,X)
  Y <- Y - Q %*% crossprod(Q,Y)
  
  bk <- backend.perm(X, Y, n.eff=n - q$rank)
  
  py <- pvals(bk, 1:dx)
  px <- pvals(bk, (1:5) + dx)
  
  pv <- c(py, px)
  expect_gte(ks.test(pv, punif)$p.value, alpha)
}

test_that("pvals with df correction", {
  set.seed(123456)
  
  test_df_correction(10, 6, 10, 1) 
  test_df_correction(40, 30, 40, 10) 
  test_df_correction(10, 6, 10, 1) 
  test_df_correction(10, 6, 10, 2) 
  test_df_correction(10, 6, 10, 3) 
})

test_that("pvals_singleton is uniform", {
  dy <- 100
  dx <- 1
  n <- 10
  alpha <- 0.01
  
  set.seed(123456)
  
  X <- matrix(rnorm(dx*n), nrow=n, ncol=dx)
  Y <- matrix(rnorm(dy*n), nrow=n, ncol=dy)
  
  bk <- backend.perm(X, Y)
        
  pv <- pvals_singleton(bk, 1)
  
  expect_equal(length(pv), dy)
  expect_gte(ks.test(pv, punif)$p.value, alpha)
})