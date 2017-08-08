pvalsCpp <- function (B) {
  
  if (min(B) > dx) {
    
    #Test X
    B <- B - dx
    cors <- if (calc_full_cor) { 
      full_xy_cor[ , B, drop = FALSE]
    } else {
      cor(X, Y[ , B])
    }
    return(pvalsC(X, Y[ , B, drop = FALSE], X4ColSum, X2, X3, cors))
    
  } else {
    
    #Test Y
    cors <- if (calc_full_cor) {
      t(full_xy_cor[B, , drop = FALSE])
    } else {
      cor(Y, X[ , B])
    }
    return(pvalsC(Y, X[ , B, drop = FALSE], Y4ColSum, Y2, Y3, cors))
    
  }
  
}

pvalsR <- function (B) {
  
  if (length(B) == 0)
    return(integer(0))
  
  test_X <- min(B) > dx
  nFixd <- length(B)
  
  if (test_X) {
    
    # Getting fixed matrix
    fixdIndx <- match(B, Yindx)
    fixdMat <- Y[ , fixdIndx, drop = FALSE]
    
    # Calculating the variances
    {
      # General calcs
      xyCors <- if (calc_full_cor) {
        full_xy_cor[ , fixdIndx, drop = FALSE]
      } else {
        cor(X, Y[ , fixdIndx, drop = FALSE])
      }
      #y4 <- colSums(X^4)
      xRowSum <- rowSums(fixdMat)
      xRowSum2 <- tcrossprod(xyCors, fixdMat^2)
      
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
    }
    
    
  } else {
    
    # Getting indices
    fixdIndx <- match(B, Xindx)
    fixdMat <- X[ , fixdIndx, drop = FALSE]
    
    # Calculating the variances
    {
      # General calcs
      xyCors <- if (calc_full_cor) { 
        full_xy_cor[fixdIndx, , drop = FALSE]
      } else {
        cor(X[,fixdIndx, drop = FALSE], Y)
      }
      xyCors <- t(xyCors)
      #y4 <- colSums(Y^4)
      xRowSum <- rowSums(fixdMat)
      xRowSum2 <- tcrossprod(xyCors, fixdMat^2)
      
      # Calc for star 1
      star1 <- crossprod(Y2, xRowSum^2)
      
      # Calc for star 2
      star2 <- Y4ColSum * rowSums(xyCors)^2
      
      # Calc for star 3
      star3 <- 2 * rowSums(xyCors) * colSums(Y2 * t(xRowSum2))
      
      # Calc for star 4
      star4 <- rowSums(xRowSum2^2)
      
      # Calc for dagger 1
      dagger1 <- rowSums(xyCors) * crossprod(Y3, xRowSum)
      
      # Calc for dagger 2
      dagger2 <- colSums(xRowSum * t(xRowSum2) * Y)
    }
    
  }
  
  allvars <- (star1 + 0.25 * (star2 + star3 + star4) - dagger1 - dagger2) / 
    (n - 1)
  corsums <- as.vector(rowSums(xyCors))
  zstats <- sqrt(n) * corsums / sqrt(allvars)
  if (twoSided) {
    pvals <- 2 * pnorm(abs(zstats), lower.tail = FALSE)
  } else {
    pvals <- pnorm(zstats, lower.tail = FALSE)
  }
  return(pvals)
  
}

pvalsR_chisq <- function (B) {
  
  if (length(B) == 0)
    return(integer(0))
  
  test_X <- min(B) > dx
  nFixd <- length(B)
  
  if (test_X) {
    
    # Getting fixed matrix
    fixdIndx <- match(B, Yindx)
    fixdMat <- Y[ , fixdIndx, drop = FALSE]
    
    if (nFixd > n) {
      
      # Do n^2|B| computation
      
    } else {
       
      # Do n|B|^2 computation
      
    }
    
  } else {
    
    # Getting indices
    fixdIndx <- match(B, Xindx)
    fixdMat <- X[ , fixdIndx, drop = FALSE]
    
    if (nFixd > n) {
      
      # Do n^2|B| computation
      
    } else {
       
      # Do n|B|^2 computation
      
    }
    
  }
  
  # Compute as and bs
  corsums <- as.vector(rowSums(xyCors))
  zstats <- sqrt(n) * corsums / sqrt(allvars)
  if (twoSided) {
    pvals <- 2 * pnorm(abs(zstats), lower.tail = FALSE)
  } else {
    pvals <- pnorm(zstats, lower.tail = FALSE)
  }
  return(pvals)
  
}