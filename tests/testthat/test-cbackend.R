source("sim_eQTL_network.R")

context("Checking C_backend")

set.seed(12344)

sim <- sim_eQTL_network(make_param_list(cmin=3, cmax=5, b=2))


dx <- ncol(sim$X) #16
dy <- ncol(sim$Y) #14
n <- nrow(sim$Y)

R <- cor(sim$X, sim$Y)

test_that("cors.base is correct when cache=0", {
  bk <- new(CorBackend, sim$X, sim$Y, 0, 0)
  expect_equal(bk$getCor(1:4 + dx), cor(sim$X, sim$Y[,1:4]))
  expect_equal(bk$getCor(4:6), cor(sim$Y, sim$X[, 4:6]))
})

test_that("cors.base is correct when cache != 0", {
  bk <- new(CorBackend, sim$X, sim$Y, 200, 200)
  
  expect_equal(bk$getCor(1:4 + dx), cor(sim$X, sim$Y[,1:4]))
  expect_equal(bk$getCor(4:6), cor(sim$Y, sim$X[, 4:6]))
  expect_equal(bk$getCor(2:5 + dx), cor(sim$X, sim$Y[,2:5]))
  expect_equal(bk$getCor(1:2), cor(sim$Y, sim$X[, 1:2]))
  expect_equal(bk$getCor(1:dx), cor(sim$Y, sim$X[, 1:dx]))
  expect_equal(bk$getCor(3:4), cor(sim$Y, sim$X[, 3:4]))
  expect_equal(bk$getCor(1:dy + dx), cor(sim$X, sim$Y[, 1:dy]))
  expect_equal(bk$getCor(3:4 + dx), cor(sim$X, sim$Y[, 3:4]))
  
  #Border cases
  expect_equal(bk$getCor(2 + dx), cor(sim$X, sim$Y[,2]))
  expect_equal(bk$getCor(2), cor(sim$Y, sim$X[, 2]))
  expect_equal(bk$getCor(dx), cor(sim$Y, sim$X[, dx]))
  expect_equal(bk$getCor(4), cor(sim$Y, sim$X[, 4]))
  expect_equal(bk$getCor(1 + dx), cor(sim$X, sim$Y[, 1]))
  expect_equal(bk$getCor(dy + dx), cor(sim$X, sim$Y[, dy]))
})