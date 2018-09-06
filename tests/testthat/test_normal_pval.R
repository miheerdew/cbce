source("sim_eQTL_network.R")

skip("Moving to cbce2")

context("Check Pvalue computation for Normal")

varcalc1_multi <- function (Ymat, Xmat) {
  
  Xmat <- scale(Xmat)
  Ymat <- scale(Ymat)
  
  n <- nrow(Xmat)
  
  if (n != nrow(Ymat))
    stop('X and Y variables must be of same dimension\n')
  
  # General calcs
  xyCors <- cor(Ymat, Xmat)
  y4 <- colSums(Ymat^4)
  xRowSum <- rowSums(Xmat)
  xRowSum2 <- tcrossprod(xyCors, Xmat^2)
  
  # Calc for star 1
  star1 <- crossprod(Ymat^2, xRowSum^2)
  
  # Calc for star 2
  star2 <- y4 * rowSums(xyCors)^2
  
  # Calc for star 3
  star3 <- 2 * rowSums(xyCors) * colSums(Ymat^2 * t(xRowSum2))
  
  # Calc for star 4
  star4 <- rowSums(xRowSum2^2)
  
  # Calc for dagger 1
  dagger1 <- rowSums(xyCors) * crossprod(Ymat^3, xRowSum)
  
  # Calc for dagger 2
  dagger2 <- colSums(xRowSum * t(xRowSum2) * Ymat)
  
  return((star1 + 0.25 * (star2 + star3 + star4) - dagger1 - dagger2) / (n - 1))
  
  
}
pval_calc <- function(X, Y){
  #Pvalue against set X
  R <- cor(X, Y)
  n <- nrow(X)
  pnorm(sqrt(n)*colSums(R)/sqrt(varcalc1_multi(Y, X)), lower.tail = FALSE)
}

sim <- sim_eQTL_network(make_param_list(cmin=3, cmax=5, b=2))
bk <- backend.normal(sim$X, sim$Y)
dx <- ncol(sim$X)
dy <- ncol(sim$Y)

test_that("Normal pval matches for testing against entire X set", {
  expect_equal(pvals(bk, 1:dx), pval_calc(sim$X, sim$Y)) 
})

test_that("Normal pval matches for testing against entire Y set", {
  expect_equal(pvals(bk, 1:dy + dx), pval_calc(sim$Y, sim$X))
})

test_that("Normal pval matches for testing against subset of X set", {
  n <- 1:(dx/2)
  expect_equal(pvals(bk, n), pval_calc(sim$X[ , n, drop=FALSE], sim$Y)) 
  n <- 5:7
  expect_equal(pvals(bk, n), pval_calc(sim$X[ , n, drop=FALSE], sim$Y)) 
  n <- 7
  expect_equal(pvals(bk, n), pval_calc(sim$X[ , n, drop=FALSE], sim$Y)) 
})

test_that("Normal pval matches for testing against subset of Y set", {
  n <- 1:(dy/2)
  expect_equal(pvals(bk, n + dx), pval_calc(sim$Y[ , n, drop=FALSE], sim$X)) 
  n <- 4:5
  expect_equal(pvals(bk, n + dx), pval_calc(sim$Y[ , n, drop=FALSE], sim$X)) 
  n <- 4
  expect_equal(pvals(bk, n + dx), pval_calc(sim$Y[ , n, drop=FALSE], sim$X)) 
})