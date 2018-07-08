context("Test cbce new")

source("sim_eQTL_network.R")
library(bmdmetrics)
library(rlist)
library(pipeR)

set.seed(1234556)

#High thresh
THRESH1 <- 0.90
#Resasonable thresh
THRESH2 <- 0.6

sim1 <- sim_eQTL_network(make_param_list(cmin=10, cmax=30, b=5, bgmult=0.05))
sim2 <- sim_eQTL_network(make_param_list(cmin=5, cmax=20, b=10, bgmult=0.05))
sim3 <- sim_eQTL_network(make_param_list(cmin=5, cmax=40, b=10, bgmult=0.1))

comms <- function(res) {
  res$extract_res %>>% 
    list.filter(success) %>>%
      list.map(StableComm)
}
report <- function(res1, sim, res2=NULL, wt.x=0.5, wt.y = 1-wt.x, weights=function(nums) nums > .9) {
  com1 <- comms(res1)
  if(is.null(res2)) {
    com2 <- comms_sim(sim)
  } else {
    com2 <- comms(res2)
  }
  compare(com1, com2, wt.x=wt.x, wt.y=wt.y, weights = weights)
}

check_sim <- function(sim, 
                      thresh1 = THRESH1, 
                      thresh2a = THRESH1,
                      thresh2b = THRESH2,
                      ...) {
  args <- list(...)
  
  res <- cbce2(sim$X, sim$Y, ...)
  
  rep <- report(res, sim)
  
  #TYPE1
  expect_true(all(rep$report1$coverage >= thresh1))
  #TYPE2
  expect_true(mean(rep$report2$coverage >= thresh2a) >= thresh2b)
}

check_results_are_almost_same <- function(res1, res2, sim, 
                                          thresh1=THRESH1,
                                          thresh2=THRESH1){
  rep <- report(res1, sim, res2)
  expect_gte(mean(rep$report1$closest_match >= thresh1), thresh2)
  expect_gte(mean(rep$report2$closest_match >= thresh1), thresh2)
}

test_that("Checking sim for cbce2", {
  check_sim(sim1)
  check_sim(sim2)
  check_sim(sim3, thresh2a=THRESH2)
})

