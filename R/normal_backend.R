#The backend for the normal test

#'@describeIn backend Constructor for one sided sum of correlations under uncorrelated Gene ans SNP sets (the weak null).
#'@inheritParams backend.base
#'@export
backend.normal <- function(X, Y, calc_full_cor=FALSE) {
  p <- backend.base(X, Y, calc_full_cor)
  p$two_sided = FALSE
  p$normal_vector_pval = function(v) stats::pnorm(sum(v), sd=sqrt(length(v)), lower.tail=FALSE)

  #p$X and p$Y are scaled version of X and Y.
  extra.precomp <- list(
                  X2 = p$X**2, X3 = p$X**3, X4ColSum = colSums(p$X**4),
                  Y2 = p$Y**2, Y3 = p$Y**3, Y4ColSum = colSums(p$Y**4)
                 )
  precomp <- append(p, extra.precomp)
  class(precomp) <- c("normal", class(p))
  precomp
}

#'@describeIn backend Constructor for two sided sum of correlations under uncorrelated Gene and SNP sets (the weak null).
#'@inheritParams backend.base
#'@export
backend.normal_two_sided <- function(X, Y, calc_full_cor=FALSE) {
  precomp <- backend.normal(X, Y, calc_full_cor)
  precomp$two_sided = TRUE
  # TODO: change precomp$normal_vector_pval
  class(precomp) <- c("normal_two_sided", class(precomp))
  precomp
}

zstats <- function(p, B) {
  #p - precomputations
  #Return the vector of z-statistics against set B.

  #We would like to pretend that B is a subset of Y nodes.
  if (min(B) > p$dx) {
    # Here that is actually the case.
    testY <- TRUE
    X2 <- p$X2
    X3 <- p$X3
    X <- p$X
    X4ColSum <- p$X4ColSum
    Y <- p$Y[, B - p$dx, drop = FALSE]
  } else {
    # Reverse the roles if that is not the case
    testY <- FALSE
    X2 <- p$Y2
    X3 <- p$Y3
    X <-  p$Y
    X4ColSum <- p$Y4ColSum
    Y <- p$X[, B, drop = FALSE]
  }

  xyCors <- cors(p, B)
  n <- p$n

  # General calcs

  #y4 <- colSums(X^4)
  xRowSum <- rowSums(Y)
  xRowSum2 <- tcrossprod(xyCors, Y^2)

  # Calc for star 1
  star1 <- crossprod(X2, xRowSum^2)

  # Calc for star 2
  star2 <- X4ColSum * rowSums(xyCors)^2

  # Calc for star 3
  star3 <- 2 * rowSums(xyCors) * colSums(X2 * t(xRowSum2))

  # Calc for star 4
  star4 <- rowSums(xRowSum2^2)

  # Calc for dagger 1
  dagger1 <- rowSums(xyCors) * crossprod(X3, xRowSum)

  # Calc for dagger 2
  dagger2 <- colSums(xRowSum * t(xRowSum2) * X)

  #Finally the variance.
  allvars <- (star1 + 0.25 * (star2 + star3 + star4) - dagger1 - dagger2) /
    (n - 1)
  corsums <- as.vector(rowSums(xyCors))
  res <- sqrt(n) * corsums / sqrt(allvars)

  if(testY) {
    res[p$maskX] <- NA
  } else {
    res[p$maskY] <- NA
  }
  res
}

#' @describeIn pvals implementation for the one sided sum of correlations under uncorrelated Gene and SNP sets (the weak null)
#' @export
pvals.normal <- function(bk, B, thresh.alpha) {
  stats::pnorm(zstats(bk, B), lower.tail = FALSE)
}

#' @describeIn pvals implementation for the two sided sum of squared correlations under uncorrelated Gene and SNP sets (the weak null)
#' @export
pvals.normal_two_sided <- function(bk, B, thresh.alpha) {
  2 * stats::pnorm(abs(zstats(bk, B)), lower.tail = FALSE)
}

#' @describeIn pvals_quick implementation for the sum of squared correlations under the strong null.
#' @export
pvals_quick.normal <- function(bk, B) {
  R <- cors(bk, B)
  if (bk$two_sided) {
    2 * stats::pnorm(abs(sqrt(bk$n) * rowSums(R)), sd=sqrt(length(B)), lower.tail = FALSE)
  } else {
    stats::pnorm(sqrt(bk$n) * rowSums(R), sd=sqrt(length(B)), lower.tail = FALSE)
  }
}
