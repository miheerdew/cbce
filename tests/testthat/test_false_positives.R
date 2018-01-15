context("Check false positives")

source("sim_eQTL_network.R")
library(bmdmetrics)

set.seed(1234556)

thresh1 <- 0.9
thresh2 <- 0.5

sim1 <- sim_eQTL_network(make_param_list(cmin=10, cmax=30, b=5, bgmult=0.05))
sim2 <- sim_eQTL_network(make_param_list(cmin=5, cmax=20, b=10, bgmult=0.05))

report <- function(res1, sim, wt.x=0.5, wt.y = 1-wt.x, weights=function(nums) nums > .9) {
  com1 <- comms(res1, sim$dx)
  com2 <- comms_sim(sim)
  compare(com1, com2, wt.x=wt.x, wt.y=wt.y, weights = weights)
}

check_sim <- function(sim, ...) {
  res <- cbce(sim$X, sim$Y, ...)
  rep <- report(res, sim)
  #TYPE1
  expect_true(all(rep$report1$coverage > thresh1))
  #TYPE2
  expect_true(mean(rep$report2$coverage > thresh1) > thresh2)
}

test_that("Checking sim for chisq", {
  check_sim(sim1, backend = 'chisq')
  check_sim(sim2, backend = 'chisq')
})

test_that("Checking sim for Normal", {
  check_sim(sim1, backend = 'normal')
  check_sim(sim2, backend = 'normal')
})