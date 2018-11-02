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
  expect_equal(fast_bh_beta(samples, alpha, a, b, lower), 
               bh_reject(pbeta(samples, a, b, lower.tail = lower), alpha))
  
  lower <- FALSE
  expect_equal(fast_bh_beta(samples, alpha, a, b, lower), 
               bh_reject(pbeta(samples, a, b, lower.tail = lower), alpha))
})

test_that("BH-reject chisq", {
  df <- 40
  alpha <- 0.5
  samples <- rchisq(n, df)
  
  lower <- TRUE

  expect_equal(fast_bh_chisq(samples, alpha, df, lower), 
               bh_reject(pchisq(samples, df, lower.tail = lower), alpha))
  
  lower <- FALSE
  expect_equal(fast_bh_chisq(samples, alpha, df, lower), 
               bh_reject(pchisq(samples, df, lower.tail = lower), alpha))
})