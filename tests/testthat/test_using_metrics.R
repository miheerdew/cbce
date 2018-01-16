context("Check false positives")

source("sim_eQTL_network.R")
library(bmdmetrics)
library(rlist)
library(pipeR)

set.seed(1234556)

THRESH1 <- 0.9
THRESH2 <- 0.6

sim1 <- sim_eQTL_network(make_param_list(cmin=10, cmax=30, b=5, bgmult=0.05))
sim2 <- sim_eQTL_network(make_param_list(cmin=5, cmax=20, b=10, bgmult=0.05))
sim3 <- sim_eQTL_network(make_param_list(cmin=5, cmax=40, b=10, bgmult=0.1))

report <- function(res1, sim, res2=NULL, wt.x=0.5, wt.y = 1-wt.x, weights=function(nums) nums > .9) {
  com1 <- comms(res1, sim$dx)
  if(is.null(res2)) {
    com2 <- comms_sim(sim)
  } else {
    com2 <- comms(res2, sim$dx)
  }
  compare(com1, com2, wt.x=wt.x, wt.y=wt.y, weights = weights)
}

check_sim <- function(sim, 
                      capture = TRUE, 
                      thresh1 = THRESH1, 
                      thresh2a = THRESH1,
                      thresh2b = THRESH2,
                      ...) {
  args <- list(...)
  
  if(capture) {
    out <- capture.output({ res <- cbce(sim$X, sim$Y, ...) })
  } else {
    res <- cbce(sim$X, sim$Y, ...)
  }
  
  rep <- report(res, sim)
  
  #TYPE1
  expect_true(all(rep$report1$coverage >= thresh1))
  #TYPE2
  expect_true(mean(rep$report2$coverage >= thresh2a) >= thresh2b)
  if("mask_extracted" %in% args && args[["mask_extracted"]]) {
    #Check that communities are disjoint.
    expect_false(res %>>%
                  list.filter(exists("StableComm") && length(StableComm) > 0) %>>%
                    unlist %>>%
                      anyDuplicated)
                  
  }
}


test_that("Checking sim for chisq", {
  check_sim(sim1, backend = 'chisq')
  check_sim(sim2, backend = 'chisq_fast')
})

test_that("Checking sim masked", {
  #TODO: Why are these test not passing?
  check_sim(sim2, thresh2a=0.85, backend = 'chisq', mask_extracted=TRUE)
  check_sim(sim3, thresh2a=0.85, backend = 'chisq_fast', mask_extracted=TRUE)
  check_sim(sim3, thresh2a=0.85, backend = 'normal', mask_extracted=TRUE)

})

test_that("Checking sim for Normal", {
  check_sim(sim1, backend = 'normal')
  check_sim(sim2, backend = 'normal')
})

test_that("chisq approximation gives the same result", {
  sim <- sim2
  out <- capture.output({ res <- cbce(sim$X, sim$Y, backend = 'chisq') })
  out <- capture.output({ res_fast <- cbce(sim$X, sim$Y, backend = 'chisq_fast') })
  rep <- report(res, sim, res_fast)
  expect_gt(min(rep$report1$closest_match), THRESH1)
  expect_gt(min(rep$report2$closest_match), THRESH1)
})