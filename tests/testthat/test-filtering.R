context("Test filtering")

source("sim_eQTL_network.R")
source("metrics.R")

library(rlist)
library(pipeR)

set.seed(1234556)

sim1 <- sim_eQTL_network(make_param_list(cmin=10, cmax=30, b=5, bgmult=0.05))

test_that("Filtering works", {
  res <- cbce(sim1$X, sim1$Y)
  expect_lte(abs(length(res$comms.fil) - length(sim1$bms)), 2)
})