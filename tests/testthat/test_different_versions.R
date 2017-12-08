source("metrics.R")
source("sim_eQTL_network.R")

sim <- sim_eQTL_network(make_param_list(cmin=5, cmax=50, b=10, bgmult=0.05))

#c1 <- cbce(sim$X, sim$Y, backend = 'chisq', init_method = "conservative-BH" )
#c2 <- cbce(sim$X, sim$Y, backend = 'chisq', init_method = "no-multiple-testing")
s1 <- cbce(sim$X, sim$Y, backend = 'normal', init_method = "conservative-BH" )
s2 <- cbce(sim$X, sim$Y, backend = 'normal', init_method = "BH-0.5" )
is1 <- cbce(sim$X, sim$Y, backend = 'indepNormal', init_method = "conservative-BH" )
ic2 <- cbce(sim$X, sim$Y, backend = 'indepChiSq', init_method = "no-multiple-testing")

report <- function(res1, res2=NULL, wt.x=0.5, wt.y = 1-wt.x, weights=function(nums) nums > .9) {
 com1 <- comms(res1, sim$dx)
 if (is.null(res2)) {
   com2 <- comms_sim(sim)
 } else {
   com2 <- comms(res2, sim$dx)
 }
 compare(com1, com2, wt.x=wt.x, wt.y=wt.y, weights = weights)
}

report_s2 <- report(s2)
report_s1 <- report(s1)
report_ic2 <- report(ic2) 
report_is1 <- report(is1) 

report_s1_ic2 <- report(s1, ic2)
