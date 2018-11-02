context("test-cutils.R")


bh_reject <- function(pvals, alpha) {
  
  #pvals may now have NAs, which should not be selected 
  m <- length(pvals)
  if (m == 0) return(integer(0))
  
  ranks <- rank(pvals, ties.method = "first")
  
  pvals_adj <-  m * pvals / ranks
  
  if (any(pvals_adj <= alpha, na.rm = TRUE)) {
    candidates <- which(pvals_adj <= alpha)
    thres <- max(pvals[candidates])
    return(which(pvals <= thres))
  } else {
    return(integer(0))
  }
  
}


set.seed(12345)
n <- 1000;

test_that("BH-reject beta", {
  a <- 30
  b <- 40
  alpha <- 0.5
  
  samples <- rbeta(n, a, b)
  
  lower <- TRUE
  expect_equal(cbce:::fast_bh_beta(samples, alpha, a, b, lower), 
               bh_reject(pbeta(samples, a, b, lower.tail = lower), alpha))
  
  lower <- FALSE
  expect_equal(cbce:::fast_bh_beta(samples, alpha, a, b, lower), 
               bh_reject(pbeta(samples, a, b, lower.tail = lower), alpha))
})

test_that("BH-reject chisq", {
  df <- 40
  alpha <- 0.5
  samples <- rchisq(n, df)
  
  lower <- TRUE

  expect_equal(cbce:::fast_bh_chisq(samples, alpha, df, lower), 
               bh_reject(pchisq(samples, df, lower.tail = lower), alpha))
  
  lower <- FALSE
  
  expect_equal(cbce:::fast_bh_chisq(samples, alpha, df, lower), 
               bh_reject(pchisq(samples, df, lower.tail = lower), alpha))
})

test_that("Jac Matrix works", {
  bimods <- list(
    list(x=1:3, y=5:3),
    list(y=5:7, x=2:5),
    list(y=1:7, x=3:20),
    list(x=2:4, y=4:7)
  )
  
  expect_identical(cbce:::jacc_matrix(bimods), cbce:::jacc_matrix_c(bimods))
})

test_that("Effective num works", {
  bimods <- list(
    list(x=1:3, y=5:3),
    list(y=5:7, x=2:5),
    list(y=1:7, x=3:20),
    list(x=2:4, y=4:7)
  )
  
  expect_identical(cbce:::effective.num1(bimods), cbce:::effective_num_c(bimods))
  
  expect_identical(cbce:::effective.num1(list()), cbce:::effective_num_c(list()))
  
})