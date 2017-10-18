source("sim_eQTL_network.R")

sim <- sim_eQTL_network(make_param_list(cmin=3, cmax=5))
res1 <- cbce(sim$X, sim$Y)
res3 <- cbce(sim$X, sim$Y, calc_full_cor = TRUE)
res2 <- cbceNW_c(sim$X, sim$Y)

all.equal(res1$report[[4]], res2$report[[4]])
list.find(1:length(res1$extract_res), !identical(res2$extract_res[[.]]$StableComm, res1$extract_res[[.]]$StableComm))
