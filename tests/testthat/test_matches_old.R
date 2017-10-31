source("sim_eQTL_network.R")
source("cbce_old.R")
library(rlist)
library(pipeR)

context("Test ChiSq")

set.seed(123345)

areIdentical <- function(res1, res2) {
 #Check extraction orders are the same
 ord1 <- res1$report[[4]]
 ord2 <- res2$report[[4]]
 
 identical(ord1, ord2) && 
 {
   s1 <- list.map(res1$extract_res[ord1], if (exists("StableComm")) StableComm else NA)
   s2 <- list.map(res2$extract_res[ord2], if (exists("StableComm")) StableComm else NA)
   all(list.zip(s1 = s1,s2 = s2) %>>% list.mapv(setequal(s1, s2)))
 }
}

test_that("ChiSq gives the same extraction", {
  
sim <- sim_eQTL_network(make_param_list(cmin=3, cmax=5, b=2))
res1 <- cbce(sim$X, sim$Y)
res2 <- cbceNW_c(sim$X, sim$Y, twoSided = TRUE)
expect_true(areIdentical(res1, res2))

sim <- sim_eQTL_network(make_param_list(cmin=3, cmax=5, b=8))
set.seed(12345)
res1 <- cbce(sim$X, sim$Y, exhaustive = TRUE)
set.seed(12345)
res2 <- cbceNW_c(sim$X, sim$Y, twoSided = TRUE, exhaustive=TRUE)
expect_true(areIdentical(res1, res2))

sim <- sim_eQTL_network(make_param_list(cmin=3, cmax=5, b=5))
res1 <- cbce(sim$X, sim$Y, start_nodes = 1:ncol(sim$X))
res2 <- cbceNW_c(sim$X, sim$Y, twoSided = TRUE, start_nodes = 1:ncol(sim$X))
expect_true(areIdentical(res1, res2))
})
