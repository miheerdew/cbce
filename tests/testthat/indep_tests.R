context("Naive backend works")
library(rlist)

dx <- 1000
dy <- 1000
n <- 100

X <- matrix(rnorm(dx*n), ncol=dx)
Y <- matrix(rnorm(dy*n), ncol=dy)

res1 <- cbce(X, Y, backend="indepChiSq", init_method = "no-multiple-testing")
res2 <- cbce(X, Y, backend="indepNormal", init_method = "no-multiple-testing", alpha = 0.05)
res2b <- cbce(X, Y, backend="indepNormal", init_method = "BH-0.5", alpha = 0.05)
res2c <- cbce(X, Y, backend="indepNormal", init_method = "non-conservative-BH", alpha = 0.05)
res3 <- cbce(X, Y, backend="normal", init_method = "no-multiple-testing", alpha = 0.05)
res4 <- cbce(X, Y, backend="chisq", init_method = "no-multiple-testing", alpha = 0.05)

bkIndepNorm <- backend.indepNormal(X, Y)
bkIndepChi <- backend.indepChiSq(X, Y)

p <- function(reps) {
  sets <- dx %/% reps
  A <- list.map(1:sets, ((. - 1)*reps + 1):(.*reps))
  pv <- function(bk) {
    unlist(list.map(A, pvals(bk, .)))
  }
  return(pv)
}

pv <- p(20)
pIndepChi <- pv(bkIndepChi)
pIndepNorm <- pv(bkIndepNorm)
pNorm <- pv(bkNorm)
pChi <- pv(bkChi)
plot <- function(ps, ps2=runif(length(ps)), ...) {
 qqplot(log(ps), log(ps2), ...)
 abline(0,1)
}

plot(pIndepNorm)
