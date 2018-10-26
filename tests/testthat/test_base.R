source("sim_eQTL_network.R")

context("Checking base_backend")

set.seed(12344)

sim <- sim_eQTL_network(make_param_list(cmin=3, cmax=5, b=2))


dx <- ncol(sim$X) #16
dy <- ncol(sim$Y) #14
n <- nrow(sim$Y)

R <- cor(sim$X, sim$Y)

test_that("cors.base is correct when cache=0", {
  bk <- backend.perm(sim$X, sim$Y, cache.size=0)
  expect_equal(cors(bk, 1:4 + dx), cor(sim$X, sim$Y[,1:4]))
  expect_equal(cors(bk, 4:6), cor(sim$Y, sim$X[, 4:6]))
})

test_that("cors.base is correct when cache != 0", {
  bk <- backend.base(sim$X, sim$Y, cache.size=200)
  
  expect_equal(cors.base(bk, 1:4 + dx), cor(sim$X, sim$Y[,1:4]))
  expect_equal(cors.base(bk, 4:6), cor(sim$Y, sim$X[, 4:6]))
  expect_equal(cors.base(bk, 2:5 + dx), cor(sim$X, sim$Y[,2:5]))
  expect_equal(cors.base(bk, 1:2), cor(sim$Y, sim$X[, 1:2]))
  expect_equal(cors.base(bk, 1:dx), cor(sim$Y, sim$X[, 1:dx]))
  expect_equal(cors.base(bk, 3:4), cor(sim$Y, sim$X[, 3:4]))
  expect_equal(cors.base(bk, 1:dy + dx), cor(sim$X, sim$Y[, 1:dy]))
  expect_equal(cors.base(bk, 3:4 + dx), cor(sim$X, sim$Y[, 3:4]))
})