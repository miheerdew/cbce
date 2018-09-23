context("Test summarizing p-val")

dx <- 10
dy <- 25
n <- 20

nsim <- 10^4

mvn <- function(n, m, rho=0.5, mu=rep(0, m)) {
  Sigma <- matrix(0, m, m)
  Sigma[1:m/2, 1:m/2] <- rho
  diag(Sigma) <- 1
  X <- mvtnorm::rmvnorm(n, sigma=Sigma, mean=mu)
}

test_that("P-value uniform under permutation", {
  X <- mvn(n, dx, mu=1:dx)
  Y <- mvn(n, dy, rho=0.1)
  
  summary <- numeric(nsim)
  
  for(i in 1:nsim) {
    summary[i] <- sum(cor(X[sample(n), ], Y)^2)
  }
  
  log.pv <- summary.pval(X, Y, summary = summary)
  testthat::expect_gt(ks.test(-log.pv, pexp)$p.value, 0.01)
})
