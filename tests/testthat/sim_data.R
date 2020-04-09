library(tibble)
library(futile.logger)


#' Efficiently generate Gaussian data that forms a single bimodule
#'
#' X matrix is generated as a gaussian with intra correlation between any two columns given by 
#' \code{rho.intra}. Then the Y matrix is obtained by linear regression of X with a 
#' \code{Ber(p)} design matrix, such that the strength of cross correlations are 
#' given by \code{rho.inter}. 
#' 
#' @param sample.size  
#' @param dim.x,dim.y The sizes of the X and Y sets of the bimodules
#' @param rho.intra,rho.inter The strength of the inter and intra correlations
#' @param d The number of X terms to select for every Y variable  \code{p*dim.x} edges.
#' 
sim_bimod <- function(
  sample.size,
  dim.x,
  dim.y,
  rho.intra,
  rho.inter,
  d
) {
  
  #----- Simulate X -------
  rho <- rho.intra
  # X ~ N(0, (1-rho) I + rho E)
  # Use X_i = \sqrt{rho} Z1  + \sqrt{1-rho} Z_i  
  Z1 <- rnorm(sample.size)
  Z <- matrix(rnorm(sample.size*dim.x), nrow=sample.size)
  X <- sqrt(1-rho)*Z + sqrt(rho)*Z1 
  
  #----- Choose a random design matrix ------
  #The number of X terms to select for every Y variable
  #d <- ceiling(p*dim.x)
  
  
  # Select d varibles from X; Just the first d for now
  # Later this vector will be randomly permuted.
  selection <- c(rep(1, d), rep(0, dim.x - d))
  
  # The design matrix
  D <- matrix(nrow=dim.y, ncol=dim.x)
  
  for(i in 1:dim.y) {
    D[i, ] <- sample(selection)
  }
  
  #----- Simulate Y using linear regression -----
  # Generate noise with appropriate variance
  c1 <- 1 - rho.intra + rho.intra*d
  
  if (c1^2/rho.inter^2 - d*c1 < 0) {
    stop("Bimodule generative constraint is violated")
  }
  
  noise.sd <- sqrt(c1^2/rho.inter^2 - d*c1)
  
  Z <- matrix(rnorm(dim.y*sample.size, sd=noise.sd), nrow=sample.size)
  
  # Simulate Y
  Y <- tcrossprod(X, D) + Z
  
  list(X=X, Y=Y)
}

#' Simulate Gaussian data containing bimodules.
#' 
#' @param sample.size The sample size
#' @param dim.x,dim.y The number of X and Y variables
#' @param num.bimod The number of bimodules to embed into the data.
#' @param background.frac The fraction of total nodes that are background nodes
#' @param edge.density The edge density for the bimodules. A vector of size num.bimod, recycled if necessary.
#' @param size.homogeniety >0, Determines how homogeneous the sizes of the bimodules are. 
#'                       Formally it represents represents the parameter of the dirichlet distribution
#'                       used to determin the bimodule sizes. If this is Inf, 
#'                       then the split will be uniform.
#' @param rho.inter,rho.intra The inter and intra correlations of the bimodule. It should be a vector of 
#'                       size num.bimod, but will be recycled if necessary.     
#' @param make_giant_component If TRUE, will add ~log(num.bimod) spurious correlations between
#'                       pairs of bimodules via extra Y variables that connect random X variables.
sim_data <- function(
  sample.size,
  rho.inter,   
  dim.x,
  dim.y=dim.x,
  rho.intra = 0,
  num.bimod = 1,
  background.frac = 0.5,
  edge.density = 0.5,
  size.homogeneity = 1,
  make_giant_component = FALSE
) {
  # --- recycle args ----
  rho.intra <- rep(rho.intra, length.out=num.bimod)
  rho.inter <- rep(rho.inter, length.out=num.bimod)
  ps <- rep(edge.density, length.out=num.bimod)
  
  
  # --- Get the number of non-background nodes ---
  m.x <- ceiling(dim.x*(1-background.frac))
  m.y <- ceiling(dim.y*(1-background.frac))
  
  # ----- Generate bimod sizes --------
  #The default distribution of sizes
  size.dist <- matrix(1/num.bimod, nrow=2, ncol=num.bimod)
  
  #If possible, use a dirichlet distribution to determine the sizes
  if(num.bimod > 1 && size.homogeneity > 0 && size.homogeneity < Inf) {
    size.dist  <- extraDistr::rdirichlet(2, size.homogeneity * rep(1, num.bimod))
  } 
  
  # Round off to get integer boundaries. Keep all positive.
  x.sizes <- pmax(floor(size.dist[1,]*m.x), 1)
  y.sizes <- pmax(floor(size.dist[2,]*m.y), 1)
  
  #Correct for rounding errors.
  #TODO: Alternatively, correct the sizes so that they add up to m.x, m.y 
  m.x <- sum(x.sizes)
  m.y <- sum(y.sizes)
  dim.x <- max(m.x, dim.x)
  dim.y <- max(m.y, dim.y)
  
  # --- Fill the bimodules into the data matix
  X <- matrix(NA, nrow=sample.size, ncol=dim.x)
  Y <- matrix(NA, nrow=sample.size, ncol=dim.y)
  
  id.x <- id.y <- 1
  
  #Table to hold meta-data about the bimodules.
  bimod.info <- tibble(rho.intra, rho.inter, ps, 
                       id.x=list(list()), 
                       id.y=list(list()))
  
  for(i in seq_len(num.bimod)) {
    
    xs <- x.sizes[i]
    ys <- y.sizes[i]
    
    re <-  rho.inter[i]
    ra <- rho.intra[i]
    d <- ceiling(ps[i]*xs)
    
    #Finish documenting bimodule meta-data.
    bimod.info$id.x[[i]] <- ids.x <- seq(id.x, id.x + xs - 1)
    bimod.info$id.y[[i]] <- ids.y <- seq(id.y, id.y + ys - 1)
    
    #In case the bimodule constraint is violated.
    if (1+ra*(d-1) < re^2 * d) {
      flog.warn("Bimodule %d inter-correlation too high; lowering it", 0)
      bimod.info$rho.inter[i] <- re <- runif(1)*sqrt((1 + ra*(d-1))/d)
    }
    
    #Finally, simulate the bimodule
    bimod <- sim_bimod(sample.size, dim.x = xs, dim.y = ys, 
                       rho.inter = re, rho.intra = ra, d = d)
    
    X[, ids.x] <- bimod$X 
    Y[, ids.y] <- bimod$Y
    
    id.x <- id.x + xs 
    id.y <- id.y + ys
    
  }
  
  print(sum(is.na(Y)))
  
  #Ensure there is a large connected component in the output bimod-bimod graph
  spurious_edges <- list()
  first_noise_id.y <- id.y
  if (make_giant_component) {
    random_edge_probability <- 2 * log(num.bimod) / num.bimod
    extra.ys <- list()
    for(i in seq_len(num.bimod)) {
      for(j in seq_len(num.bimod)) {
        if (runif(1) < random_edge_probability) {
          #Choose random connection points
          x1 <- sample(bimod.info$id.x[[i]], 1)
          x2 <- sample(bimod.info$id.x[[j]], 1)
          
          #Add a new variable that connects x1 and x2
          rho_new <- min(mean(bimod.info$rho.inter[c(i, j)]), 1 / sqrt(2))
          std_new <- sqrt(1 / rho_new^2 - 2)
          extra.ys[[length(extra.ys) + 1]]  <-
            X[, x1] + X[, x2] + rnorm(sample.size) * std_new
          spurious_edges[[length(extra.ys)]] <-
            data.frame(
              x=c(x1, x2),
              y=c(id.y, id.y),
              bimod=c(i, j),
              stringsAsFactors=FALSE)
          id.y <- id.y + 1
        }
      }
    }
    #Add extra ys
    deficit <- id.y - dim.y - 1
    if (deficit > 0) {
      Y <- cbind(Y, matrix(NA, nrow=sample.size, ncol=deficit))
      dim.y <- id.y - 1
    }
    Y[, seq(first_noise_id.y, id.y - 1)] <- do.call(cbind, extra.ys)
  }
  bridge_edgelist <- do.call(rbind, spurious_edges)
  
  # If there is space, fill with iid white noise.
  if(id.x <= dim.x) {
    X[, seq(id.x, dim.x)] <- rnorm(sample.size*(dim.x - id.x + 1))
  }
  print(id.y)
  print(dim.y)
  print(dim.y - id.y + 1)
  if(id.y <= dim.y) {
    Y[, seq(id.y, dim.y)] <- rnorm(sample.size*(dim.y - id.y + 1))
  }
  print(sum(is.na(Y)))
  list(X=X, Y=Y, bimod.info=bimod.info, bridge_edgelist=bridge_edgelist)
}


if(FALSE) {
  data <- sim_data(100,
                   num.bimod = 3,
                   dim.x=100, dim.y = 100,
                   size.homogeneity = 10,
                   rho.inter=c(0.2, 0.3, 0.4),
                   rho.intra = c(0.1, 0.1, 0.1), 
                   edge.density = c(0.9, .3, .1))
  
  cormat <- cor(cbind(data$X, data$Y))
  
  library(ggplot2)
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
  
  data <- sim_data(100, 0.4, dim.x = 200, dim.y = 500, rho.intra=0.3, num.bimod = 3)
  
  cormat <- cor(cbind(data$X, data$Y))
  
  ggcor(cormat, fisher = FALSE, fn = "default_mat2.png", theme = cortheme,
        title = "Default Network Model", width = 11.5, height = 10)
}


