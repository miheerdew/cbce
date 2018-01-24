source("sim_eQTL_network.R")
source("cbce_old.R")
library(rlist)
library(pipeR)

context("Test ChiSq with old version")

set.seed(123345)

extraction_ind <- function(exres) {
  exres %>>% 
    list.which(!(report %in% c('dud', 'indx_clustered')))
}

check <- function(res1, res2, 
                  expect_differing_indices=FALSE) {
 #Check extraction orders are the same
 ord1 <- res1$report[[4]]
 ord2 <- res2$report[[4]]
 
 expect_identical(ord1, ord2)
 
  ex1 <- res1$extract_res[ord1]
  ex2 <- res2$extract_res[ord2]
   
  ind1 <- extraction_ind(ex1)
  ind2 <- extraction_ind(ex2)

  if(expect_differing_indices) {
    expect_false(identical(ind1, ind2))
    if(all(ind2 %in% ind1)){
      ind = ind2
    # Only allow one way inclusion.  
    #} else if (all(ind2 %in% ind1)) {
    #ind = ind2
    } else {
     #should not be reached.
     expect_true(FALSE)
    }
  } else {
    expect_identical(ind1, ind2)
    ind = ind1 
  }
  
  list.iter(ind, 
        expect_setequal(ex1[[.]]$StableComm, ex2[[.]]$StableComm)
    )
}

test_that("ChiSq gives the same extraction", {
out <- 
capture.output({  
  sim <- sim_eQTL_network(make_param_list(cmin=10, cmax=30, b=5, bgmult=0.05))
  res1 <- cbce(sim$X, sim$Y)
  res2 <- cbceNW_c(sim$X, sim$Y, twoSided = TRUE)
  check(res1, res2)
  
  sim <- sim_eQTL_network(make_param_list(cmin=10, cmax=20, b=5, bgmult=0.01))
  #Since we fixed the threshold bug...
  # The first half fails, but second half passes.
  #set.seed(12345)
  #res1a <- cbce(sim$X, sim$Y, exhaustive = TRUE, OL_thres = 0.7)
  set.seed(12345)
  res1b <- cbce(sim$X, sim$Y, exhaustive = TRUE, OL_thres = 0.5)
  set.seed(12345)
  res2 <- cbceNW_c(sim$X, sim$Y, twoSided = TRUE, exhaustive=TRUE)
  #expect_false(areIdentical(res1a, res2, sim$dx))
  
  #Indices differ because new dud scheme is more relaxed.
  check(res1b, res2, expect_differing_indices = TRUE)
  
  
  sim <- sim_eQTL_network(make_param_list(cmin=10, cmax=30, b=5, bgmult=0.05))
  res1 <- cbce(sim$X, sim$Y, start_nodes = 1:sim$dx)
  res2 <- cbceNW_c(sim$X, sim$Y, twoSided = TRUE, start_nodes = 1:sim$dx)
  check(res1, res2)
})
})
