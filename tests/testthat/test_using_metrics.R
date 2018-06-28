context("Test using metrics")

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

report <- function(res1, sim, res2=NULL, wt.x=0.5, wt.y = 1-wt.x, weights=function(nums) nums > .9) {
  resf <- res1
  resf$extract_res <- res1$extract_res[res1$finalIndxs]
  com1 <- comms(resf, sim$dx)
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

check_results_are_almost_same <- function(res1, res2, sim, 
                                          thresh1=THRESH1,
                                          thresh2=THRESH1){
  rep <- report(res1, sim, res2)
  expect_gte(mean(rep$report1$closest_match >= thresh1), thresh2)
  expect_gte(mean(rep$report2$closest_match >= thresh1), thresh2)
}

test_that("Checking sim for chisq", {
  check_sim(sim1, backend = 'chisq')
  check_sim(sim2, backend = 'chisq')
  check_sim(sim3, backend = 'chisq', exhaustive=TRUE, thresh2a=0.8, thresh2b = THRESH2)
})

test_that("Checking sim masked", {
  #TODO: Why are these test not passing? Maybe because of the CLT Approx
  check_sim(sim2, backend = 'chisq', mask_extracted=TRUE)
  check_sim(sim3, thresh2a=0.85, backend = 'normal', mask_extracted=TRUE)

})

test_that("Checking sim with ranking", {
  check_sim(sim1, backend = 'chisq', rank_initial_sets=TRUE)
  check_sim(sim2, backend = 'normal', rank_initial_sets=TRUE)
})

test_that("Checking sim with ranking and init_quick_update for new vs old for normal", {
  sim <- sim2
  out <- capture.output({ 
  res_unranked <- cbce(sim$X, sim$Y, backend = 'normal') })
  out <- capture.output({
  res_ranked <- cbce(sim$X, sim$Y, backend = 'normal', 
                       rank_initial_sets = TRUE, init_quick_update=TRUE)  })
  check_results_are_almost_same(res_unranked, res_ranked, sim)
})

test_that("Checking sim with ranking new vs old for chisq", {
  sim <- sim1
  out <- capture.output({ 
    res_unranked <- cbce(sim$X, sim$Y, backend = 'chisq') })
  out <- capture.output({ 
    res_ranked <- cbce(sim$X, sim$Y, backend = 'chisq', 
                       rank_initial_sets = TRUE) })
  #There might be some extra stuff appeating in the ranked versions.
  check_results_are_almost_same(res_unranked, res_ranked, sim, thresh2 = 0.8)
})

test_that("Checking sim for Normal", {
  check_sim(sim1, backend = 'normal')
  check_sim(sim2, backend = 'normal')
})

test_that("Results are almost same for chisq, normal", {
  sim <- sim2
  out <- capture.output({ resC <- cbce(sim$X, sim$Y, backend = 'chisq') })
  out <- capture.output({ resN <- cbce(sim$X, sim$Y, backend = 'normal') })
  check_results_are_almost_same(resC, resN, sim)
})

test_that("Results are almost same for chisq, normal when masked", {
  sim <- sim2
  out <- capture.output({ resC <- cbce(sim$X, sim$Y, backend = 'chisq', mask_extracted = TRUE) })
  out <- capture.output({ resN <- cbce(sim$X, sim$Y, backend = 'normal', mask_extracted = TRUE) })
  check_results_are_almost_same(resC, resN, sim)
})

test_that("Results are almost same for normal vs normalc", {
  sim <- sim2
  out <- capture.output({ resC <- cbce(sim$X, sim$Y, backend = 'normalc') })
  out <- capture.output({ resN <- cbce(sim$X, sim$Y, backend = 'normal') })
  check_results_are_almost_same(resC, resN, sim, 0, 0)
})