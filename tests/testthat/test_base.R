source("sim_eQTL_network.R")

context("Checking base_backend")

set.seed(12344)
sim <- sim_eQTL_network(make_param_list(cmin=3, cmax=5, b=2))
bk <- backend.normal(sim$X, sim$Y)
bk2 <- backend.chisq(sim$X, sim$Y)


dx <- ncol(sim$X) #16
dy <- ncol(sim$Y) #14
n <- nrow(sim$Y)

R <- cor(sim$X, sim$Y)

X <- matrix(rnorm(dx*n), ncol=dx)
Y <- matrix(rnorm(dy*n), ncol=dy)
bk3 <- backend.normal(X, Y)
bk4 <- backend.chisq(X, Y)

test_that("cors.base gives correlations", {
  expect_equal(cors.base(bk, 1:8, 1:4 + dx), cor(sim$Y[,1:4], sim$X[,1:8]))
  expect_equal(cors.base(bk, 8:1+dx, 4:6), cor(sim$X[,4:6], sim$Y[,8:1]))
  expect_equal(cors.base(bk, 8:1+dx), cor(sim$X, sim$Y[,8:1]))
})

test_that("score_quick calculates sensible score", {
  expect_gt(score_quick(bk, 1:dx, 1:dy+dx), 10)
  expect_gt(score_quick(bk2, 1:dx, 1:dy+dx), 50)
  expect_equal(score_quick(bk, 1:dx, 1:dy+dx), score_quick(bk, 1:dy + dx, 1:dx))
  expect_equal(score_quick(bk2, 1:dx, 1:dy+dx), score_quick(bk2, 1:dy + dx, 1:dx))
})

test_that("score_quick low under null", {
  expect_lt(score_quick(bk3, 1:dx, 1:dy+dx), 2)
  expect_lt(score_quick(bk4, 1:dx, 1:dy+dx), 2)
})
