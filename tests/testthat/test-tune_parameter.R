source("sim_eQTL_network.R")

set.seed(1234556)

sim1 <- sim_eQTL_network(make_param_list(cmin=10, cmax=30, b=5, bgmult=0.05))


test_that("Parameter tuning works", {
  capture.output({
  fdrs <- half_permutation_fdr(sim1$X, sim1$Y, 
                               alphas=c(0.01, 0.02, 0.03),
                               num.sims=5)
  })
  expect_lt(max(fdrs), 0.05)
})

test_that("Parameter timing out", {
  capture.output({
    fdrs <- half_permutation_fdr(sim1$X, sim1$Y, 
                               alphas=c(0.01),
                               num.sims=2,
                               method=function(...) Sys.sleep(10),
                               timeout=5)
  })
  expect_true(all(is.na(fdrs)))
})

test_that("Parameter tuning with timeout", {
  capture.output({
    fdrs <- half_permutation_fdr(sim1$X, sim1$Y, 
                                 alphas=c(0.01, 0.02, 0.03),
                                 num.sims=5,
                                 timeout=5)
  })
  expect_lt(max(fdrs), 0.05)
})
