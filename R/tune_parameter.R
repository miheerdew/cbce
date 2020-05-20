#' Estimate the false discovery rate for the method across a 
#' variety of parameters.
#' 
#' For each parameter value, generate \code{num.sims} half-permuted 
#' datasets, and estimate the half-permutation FDR.
#' 
#' @param alphas List of parameters to run the method on.
#' @param X,Y,cov The original data that we intend to run the method 
#'                on.
#' @param num.sims The number of simulation to evaluate over.
#' @param method The method to run. The method should take 
#'               arguments \code{X, Y, cov, alpha} and return 
#'               a collection of bi-modules. 
#' @param fdr The type of false discovery metric to use. Can be one of 'all.pairs' or 'imp.pairs'.
#'
#' @param timeout The maximum amount of time in seconds to let each method run, 
#' before stopping it. If set to Inf or NULL, the method will never be timed out.
#' 
#' @return 
#' Returns a numeric vector of size length(alphas) that have the fdrs.
#'@keywords internal
#'@export
half_permutation_fdr <- function(X, Y, alphas,  
                       num.sims, method=cbce.fast.fil, cov=NULL,
                       fdr='imp.pairs', filter=FALSE, timeout=Inf) {
  
  if(is.null(timeout)) timeout <- Inf
  
  if(timeout < Inf) {
    if (!requireNamespace("R.utils", quietly = TRUE)) {
      stop("Package \"R.utils\" needed for timeout functionality to work. Please install it.",
           call. = FALSE)
    }
  }

  fdr.mat <- matrix(numeric(num.sims*length(alphas)),
                 nrow=num.sims)

  nsample <- nrow(X)
  
  for(i in seq_len(num.sims)) {
    scr.cols=list(x=NULL, y=NULL)
    
    scr.cols$y <- which(rbinom(ncol(Y), 1, 0.5) == 1)
    Y.scr <- Y
    Y.scr[, scr.cols$y] <- Y.scr[sample(nsample), scr.cols$y]

    scr.cols$x <- which(rbinom(ncol(X), 1, 0.5) == 1)
    X.scr <- X
    X.scr[, scr.cols$x] <- X.scr[sample(nsample), scr.cols$x]
    
    for(j in seq_along(alphas)) {
      print(sprintf("Data %d, alpha=%.2E", i, alphas[j]))
      if(timeout < Inf) {
        bimods <- NULL
        e <- try(
          bimods <- R.utils::withTimeout(method(X.scr, Y.scr, alphas[j], cov), 
                                        timeout=timeout, onTimeout = "error"),
          silent=FALSE,
          outFile = stdout()
        )
      } else {
        e <- NULL
        bimods <- method(X.scr, Y.scr, alphas[j], cov)
      }
      bimods <- filter_bimodules(bimods)
      fds <- purrr::map_dbl(bimods, ~
                       switch(fdr,
                              all.pairs=P.pairs(scr.cols, .),
                              imp.pairs=FDR.imp_pairs(
                                cor(
                                  X.scr[, .$x, drop=FALSE], 
                                  Y.scr[, .$y, drop=FALSE]), 
                                na.omit(match(scr.cols$x, .$x)),
                                na.omit(match(scr.cols$y, .$y))
                              )))
      
      if (inherits(e, "try-error")) {
        fdr.mat[i,j] <- NA
      } else {
        fdr.mat[i, j] <- if(length(fds) > 0) mean(fds) else 0
      }
    }
  }
  colMeans(fdr.mat)
}

#P(A|B)
P <- function(A, B) {
  length(intersect(A, B))/length(B)
}

#P(Ap|Bp) for pairs
P.pairs <- function(Ap, Bp) {
  cy <- length(intersect(Ap$y, Bp$y))
  cx <- length(intersect(Ap$x, Bp$x))
  
  sx <- length(Bp$x)
  sy <- length(Bp$y)
  
  (sx*cy + cx*sy - cx*cy)/(sx*sy)
}

FDR.imp_pairs <- function(R, scr.x, scr.y) {
  th <- connectivity_threshold(R)
  M <- abs(R) > th
  tot <- sum(M)
  fd <- sum(M[scr.x,]) + sum(M[, scr.y]) - sum(M[scr.x, scr.y])
  fd/tot
}

cbce.fast <- function(X, Y, alpha, cov=NULL) {
  cbce(X, Y, alpha, cov=cov, heuristic_search=TRUE,
       interaction=interaction_cli, start_frac = 0.5,
       filter_low_score=FALSE)$comms
}

cbce.fast.fil <- function(X, Y, alpha, cov=NULL) {
  cbce(X, Y, alpha, cov=cov, heuristic_search=TRUE,
       interaction=interaction_cli,
       filter_low_score=FALSE)$comms.fil
}

cbce.fil <- function(X, Y, alpha, cov=NULL) {
  cbce(X, Y, alpha, cov=cov, heuristic_search=TRUE,
       interaction=interaction_cli,
       filter_low_score=FALSE)$comms.fil
}

cbce.bare <- function(X, Y, alpha, cov=NULL) {
  cbce(X, Y, alpha, cov=cov,
       interaction=interaction_cli,
       filter_low_score=FALSE)$comms
}