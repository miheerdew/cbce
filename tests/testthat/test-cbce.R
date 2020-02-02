context("Test cbce")

source("sim_eQTL_network.R")
source("metrics.R")

library(rlist)
library(pipeR)

set.seed(1234556)

#High thresh
THRESH1 <- 0.90
THRESH75 <- 0.75
#Resasonable thresh
THRESH2 <- 0.6

sim1 <- sim_eQTL_network(make_param_list(cmin=10, cmax=30, b=5, bgmult=0.05))
sim2 <- sim_eQTL_network(make_param_list(cmin=5, cmax=20, b=10, bgmult=0.05))
sim3 <- sim_eQTL_network(make_param_list(cmin=5, cmax=40, b=10, bgmult=0.1))

report <- function(com1, sim, res2=NULL, wt.x=0.5, wt.y = 1-wt.x, weights=function(nums) nums > .9) {
  if(is.null(res2)) {
    com2 <- comms_sim(sim)
  } else {
    com2 <- comms(res2)
  }
  compare(com1, com2, wt.x=wt.x, wt.y=wt.y, weights = weights)
}

check_sim <- function(sim, 
                      add.cov = NULL,
                      thresh1 = THRESH1, 
                      thresh2a = THRESH1,
                      thresh2b = THRESH2,
                      thresh3 = THRESH75,
                      ...) {
  args <- list(...)
  
  if(!is.null(add.cov)) {
    sample.size <- nrow(sim$X)
    cvrt <- matrix(rnorm(add.cov*sample.size), nrow=sample.size)
    X <- sim$X + rowSums(cvrt)
    Y <- sim$Y + rowSums(cvrt)
    res <- cbce(X, Y, cov=cvrt, ...)
  } else {
    res <- cbce(sim$X, sim$Y, ...)
  }
  
  rep <- report(res$comms.fil, sim)
  
  #TYPE1
  expect_true(all(rep$report1$coverage >= thresh1))
  #TYPE2
  expect_true(mean(rep$report2$coverage >= thresh2a) >= thresh2b)
  
  n <- length(sim$bms)
  expect_lte(abs(res$summary$eff.num - n)/n, 1 - thresh3)
}

check_results_are_almost_same <- function(res1, res2, sim, 
                                          thresh1=THRESH1,
                                          thresh2=THRESH1){
  rep <- report(res1, sim, res2)
  expect_gte(mean(rep$report1$closest_match >= thresh1), thresh2)
  expect_gte(mean(rep$report2$closest_match >= thresh1), thresh2)
}

test_that("Checking sim for cbce", {
  check_sim(sim1)
  check_sim(sim2)
  check_sim(sim3, thresh2a=THRESH2,
            alpha=0.01)
})

test_that("Check heuristic_search for type1 error", {
  check_sim(sim1, heuristic_search=TRUE)
  check_sim(sim2, heuristic_search=TRUE)
  check_sim(sim3, heuristic_search=TRUE, thresh2a=THRESH2)
})

test_that("Test for correction of covariate", {
  check_sim(sim1, heuristic_search=TRUE, add.cov=3)
  check_sim(sim2, heuristic_search=TRUE, add.cov=7)
  check_sim(sim3, heuristic_search=TRUE, thresh2a=THRESH2, add.cov=10)
})

test_that("Parameter tuning works", {
  fdrs <- half_permutation_fdr(sim1$X, sim1$Y, 
                               alphas=c(0.01, 0.02, 0.03),
                               num.sims=5)
  expect_lt(max(fdrs), 0.05)
})