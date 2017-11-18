context("Checking independent backend")

dx <- 100
dy <- 150
n <- 500
alpha <- 0.01

set.seed(12345)
X <- matrix(rnorm(dx*n), ncol=dx)
Y <- matrix(rnorm(dy*n), ncol=dy)
bkNorm <- backend.indepNormal(X, Y)
bkChi <- backend.indepChiSq(X, Y)

mypvals <- function(A) {
  if(min(A) > dx){
    A <- A - dx
    R <- cor(X, Y[, A])
  } else {
    R <- cor(Y, X[,A])
  }
  k <- length(A)
  t <- sqrt(n)*rowSums(R)
  t2 <- n*rowSums(R^2)
  list(norm=pnorm(t, sd=sqrt(k), lower.tail = FALSE),
       chisq=pchisq(t2, k, lower.tail = FALSE))
}
test_equal <- function(A) {
  p1 <- mypvals(A)
  expect_equal(p1$norm, pvals(bkNorm, A))
  expect_equal(p1$chisq, pvals(bkChi, A))
}
test_uniform <- function(A){
  expect_gte(ks.test(pvals(bkNorm, A), "punif")$p.value, alpha)
  expect_gte(ks.test(pvals(bkNorm, A), "punif")$p.value, alpha)
}

test_that("Pvalues are calculated as expected", {
  test_equal(2:35)
  test_equal(101:150)
})

test_that("Pvalues are uniform", {
  test_uniform(1:100)
  test_uniform(120:150)
})

