# Function to make eQTL-network-like data
# Argument descriptions:
# n = number of samples
# b = number of bimodules
# cmin = minimum size of bimodule half
# cmax = maximum size of bimodule half
# bgmult = number of background variables / (b * (cmax - cmin) / 2)
# betamean = mean of exponential distribution for betas
# rho = intracorrelations of X variables
# p = average edge density of bimodules
# s2 = noise variance scaling
# bgb = number of background blocks to make

### naive implementation in R
mvrnormR <- function(n, mu, sigma) {
  ncols <- ncol(sigma)
  mu <- rep(mu, each = n) ## not obliged to use a matrix (recycling)
  mu + matrix(rnorm(n * ncols), ncol = ncols) %*% chol(sigma)
}
sim_eQTL_network <- function (par_list, randomizeBeta = TRUE) {
  
  
  n = par_list$n
  b = par_list$b
  cmin = par_list$cmin
  cmax = par_list$cmax
  bgmult = par_list$bgmult
  betamean = par_list$betamean
  p = par_list$p
  rho = par_list$rho
  s2 = par_list$s2
  bgb = par_list$bgb
  corNoiseX = par_list$corNoiseX
  corNoiseY = par_list$corNoiseY
  maxNoiseCor = par_list$maxNoiseCor
  
  avgsize <- cmin + (cmax - cmin) / 2
  
  # Setting noise variance
  mavg <- p * avgsize
  nv <- s2 * mavg * (1 - rho + rho * mavg) * betamean
  
  # Getting bimodule sizes
  if (cmin < cmax) {
    Xsizes <- sample(cmin:cmax, b, replace = TRUE)
    Ysizes <- sample(cmin:cmax, b, replace = TRUE)
  } else {
    Xsizes <- Ysizes <- rep(cmin, b)
  }
  
  # Getting number of background blocks
  dbg <- round(bgmult * b * avgsize)
  if (is.null(bgb)) {
    bgb <- round(dbg / avgsize)
    bgb <- max(bgb, 1)
  }
  dbg <- max(dbg, 3 * bgb)
  dx <- sum(Xsizes) + dbg; dy <- sum(Ysizes) + dbg
  
  # Initializing loop
  X <- matrix(numeric(dx * n), nrow = n); Y <- matrix(numeric(dy * n), nrow = n)
  Xindx <- 1; Yindx <- 1
  bms <- rep(list(NULL), b)
  for (i in 1:b) {
    
    # Making indices
    Xindxs <- Xindx:(Xindx + Xsizes[i] - 1)
    Yindxs <- Yindx:(Yindx + Ysizes[i] - 1)
    bms[[i]] <- c(Xindxs, dx + Yindxs)
    MindxX <- 1:Xsizes[i]
    MindxY <- (Xsizes[i] + 1):(Xsizes[i] + Ysizes[i])
    
    # Generating eQTLs
    eQTLs <- matrix(rbinom(Xsizes[i] * Ysizes[i], 1, p), ncol = Ysizes[i])
    if (randomizeBeta)
      eQTLs[eQTLs == 1] <- rexp(sum(eQTLs), rate = 1 / betamean)
    
    # Making data
    SigmaX <- (1 - rho) * diag(Xsizes[i]) + 
      matrix(rho, ncol = Xsizes[i], nrow = Xsizes[i])
    Xi <- mvrnormR(n, rep(0, Xsizes[i]), SigmaX)
    noisevec <- rnorm(n * Ysizes[i], sd = sqrt(nv))
    noisemat <- matrix(noisevec, nrow = n)
    Yi <- Xi %*% eQTLs + noisemat
    
    # Putting data in, re-setting loop
    X[ , Xindxs] <- Xi; Y[ , Yindxs] <- Yi
    Xindx <- Xindx + Xsizes[i]; Yindx <- Yindx + Ysizes[i]
        
  }
  
  if (dbg > 0) {
    
    # Making background blocks----------------------------------------------------
    
    if (bgb == 1) {
      Xbgindxs <- Ybgindxs <- list(1:dbg)
    } else {
      breakpointsX <- sort(sample(2:(dbg - 1), bgb - 1, replace = FALSE))
      breakpointsY <- sort(sample(2:(dbg - 1), bgb - 1, replace = FALSE))
    }
    if (bgb == 2) {
      
      Xbgindxs <- c(list(1:breakpointsX[1]), 
                    list(breakpointsX[bgb - 1]:dbg))
      Ybgindxs <- c(list(1:breakpointsY[1]), 
                    list(breakpointsY[bgb - 1]:dbg))
    }
    if (bgb > 2) {
      Xbgindxs <- lapply(2:(bgb - 1), 
                         function (i) breakpointsX[i - 1]:(breakpointsX[i] - 1))
      Xbgindxs <- c(list(1:breakpointsX[1]), 
                    Xbgindxs, 
                    list(breakpointsX[bgb - 1]:dbg))
      
      Ybgindxs <- lapply(2:(bgb - 1), 
                         function (i) breakpointsY[i - 1]:(breakpointsY[i] - 1))
      Ybgindxs <- c(list(1:breakpointsY[1]), 
                    Ybgindxs, 
                    list(breakpointsY[bgb - 1]:dbg))
    }
    
    # Getting intracorrs
    Xrhos <- runif(bgb, 0, maxNoiseCor)
    Yrhos <- runif(bgb, 0, maxNoiseCor)
    if (!corNoiseX) Xrhos <- rep(0, bgb)
    if (!corNoiseY) Yrhos <- rep(0, bgb)
    
    Xnoise <- Ynoise <- matrix(0, ncol = dbg, nrow = n)
    for (j in 1:bgb) {
      Xbgbj <- length(Xbgindxs[[j]]); Ybgbj <- length(Ybgindxs[[j]])
      SigXj <- matrix(rep(Xrhos[j], Xbgbj^2), ncol = Xbgbj) + 
        diag(Xbgbj) * (1 - Xrhos[j])
      SigYj <- matrix(rep(Yrhos[j], Ybgbj^2), ncol = Ybgbj) +
        diag(Ybgbj) * (1 - Yrhos[j])
      Xnoise[ , Xbgindxs[[j]]] <- mvrnormR(n, rep(0, Xbgbj), SigXj)
      Ynoise[ , Ybgindxs[[j]]] <- mvrnormR(n, rep(0, Ybgbj), SigYj)
    }
  
    X[ , Xindx:dx] <- Xnoise
    Y[ , Yindx:dy] <- Ynoise
    
  }
  
  return(list("X" = X, "Y" = Y, "bms" = bms, "dx"=ncol(X)))
  
}

make_param_list <- function (n = 200,
                             b = 10,
                             cmin = 50,
                             cmax = 100,
                             bgmult = 1,
                             betamean = 1,
                             p = 0.5,
                             rho = 0.3,
                             s2 = 4,
                             bgb = NULL,
                             corNoiseX = FALSE,
                             corNoiseY = FALSE,
                             maxNoiseCor = rho) {
  return(list("n" = n,
              "b" = b,
              "cmin" = cmin,
              "cmax" = cmax,
              "bgmult" = bgmult,
              "betamean" = betamean,
              "p" = p,
              "rho" = rho,
              "s2" = s2,
              "bgb" = bgb,
              "corNoiseX" = corNoiseX,
              "corNoiseY" = corNoiseY,
              "maxNoiseCor" = maxNoiseCor))
}


if (FALSE) {
  
  par_list <- make_param_list()
  default_sim <- sim_eQTL_network(par_list)
  cormat <- cor(cbind(default_sim$X, default_sim$Y))
  cortheme <- theme(axis.title.x = element_text(size = 20),
                    axis.title.y = element_text(size = 20),
                    axis.text.x = element_text(size = 15),
                    axis.text.y = element_text(size = 15),
                    title = element_text(size = 25),
                    legend.title = element_text(size = 20),
                    legend.text = element_text(size = 15))
  
  source("ggcor.R")
  ggcor(cormat, fisher = FALSE, fn = "default_mat.png", theme = cortheme,
        title = "Default Network Model", width = 11.5, height = 10)
}
    
    