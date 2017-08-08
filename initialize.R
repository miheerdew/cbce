initializeR <- function(cors, alpha, conserv = TRUE) {
  fischer_tranformed_cor <- atanh(cors) * sqrt(n - 3)
  if (twoSided) {
    pvals <- 2 * pnorm(abs(fischer_tranformed_cor), lower.tail = FALSE)
  } else {
    pvals <- pnorm(fischer_tranformed_cor, lower.tail = FALSE)
  }
  successes <- bh_rejectR(pvals, alpha)
  return(successes)
}

initialize1 <- function (u, Cpp = TRUE) {
  
  if (u <= dx) {
    
    # Test X
    cors <- if (calc_full_cor) { 
      t(full_xy_cor[u, ])
    } else {
      cor(Y, X[ , u])
    }
    
    if (Cpp) {
      return(initializeC(n, cors, alpha, conserv = TRUE) + dx)
    } else {
      return(initializeR(cors, alpha, conserv = TRUE) + dx)
    }
    
    
  } else {
    
    # Test Y
    cors <- if (calc_full_cor) {
      full_xy_cor[ , u - dx]
    } else {
      cor(X, Y[ , u - dx])
    }
    
    if (Cpp) {
      return(initializeC(n, cors, alpha, conserv = TRUE))
    } else {
      return(initializeR(cors, alpha, conserv = TRUE))
    }

  }
  
  return(successes)
  
}