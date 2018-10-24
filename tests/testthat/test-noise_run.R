context("test-noise_run")

set.seed(12345)

n <- 100
dx <- 50
dy <- 70


X <- matrix(rnorm(n*dx), ncol=dx)
Y <- matrix(rnorm(n*dy), ncol=dy)

test_that("No bimodule on noise", {
  res <- cbce2(X, Y)
  expect_null(res$filtered_result.df)
})

stop_after_five <- function(event, env) {
  if(event == "Main:NextExtraction") {
    env$count <- get0("count", envir=env, ifnotfound = 1) + 1
    if(env$count > 5){
      #return("browse")
      return("stop")
    }
  }
}

test_that("Run stopped", {
  res <- cbce2(X, Y, interaction=stop_after_five)
  expect_null(res$filtered_result.df)
})
