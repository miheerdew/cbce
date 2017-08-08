library(Rcpp)
library(RcppParallel)
library(foreach)
library(doParallel)
library(bmdupdate)

#sourceCpp("correlation.cpp")

cbceNW_c <- function (X, Y, alpha = 0.05, OL_thres = 0.9, tag = NULL, Cpp = FALSE, verbose = TRUE, generalOutput = TRUE,
                      updateOutput = TRUE, exhaustive = FALSE, OL_tol = Inf, Dud_tol = Inf, time_limit = 18000,
                      updateMethod = 1, inv.length = TRUE, add_rate = 1, start_nodes = NULL, pval_parallel = FALSE,
                      calc_full_cor=FALSE, loop_limit = Inf, parallel = FALSE, twoSided = FALSE) {
  
  if (FALSE) {
    alpha = 0.05
    OL_thres = 0.9
    tag = NULL
    Cpp = FALSE
    verbose = TRUE
    generalOutput = TRUE
    OL_tol = Inf
    Dud_tol = Inf
    time_limit = Inf
    updateMethod = 1
    updateOutput = TRUE
    exhaustive = FALSE
    inv.length = TRUE
    loop_limit = Inf
    add_rate = 1
    calc_full_cor = TRUE
    parallel = FALSE
    start_nodes = NULL
    pval_parallel = FALSE
    twoSided = FALSE
  }
  
  if(Cpp) {
    library(bmdCpp)
  }
  #-----------------------------------------------------------------------------
  # Setup 
  
  start_second <- proc.time()[3]
  td <- tempdir()
  cat("#-------------------\n")
  
  if (generalOutput)
    cat("Setup & pre-calculations\n")
  
  X <- scale(X); Y <- scale(Y)
  dx <- ncol(X)
  dy <- ncol(Y)
  n  <- nrow(X)
  Xindx <- 1:dx
  Yindx <- (dx + 1):(dx + dy)
  
  if(calc_full_cor){
    if (generalOutput)
      cat("Calculating full cross correlation matrix\n")
    full_xy_cor <- cor(X, Y)
  }
  bmd_obj <- new(BmdUpdater, X, Y)
  
  
  #-------------------------------------------------------------------------------
  # Auxiliary functions
  symdiff <- function (s1, s2) {
    return(union(setdiff(s1, s2), setdiff(s2, s1)))
  }
  
  jaccard <- function (s1, s2) {
    return(length(symdiff(s1, s2)) / length(union(s1, s2)))
  }
  
  filter_overlap <- function (comms, tau, inv.length = FALSE) {
    
    K <- length(comms)
    if (inv.length) {
      scores <- 1 / unlist(lapply(comms, length))
    } else {
      scores <- unlist(lapply(comms, length))
    }
    
    
    jaccard_mat0 <- matrix(0, K, K)
    for (i in 1:K) {
      for (j in 1:K) {
        jaccard_mat0[i, j] <- length(intersect(comms[[i]], comms[[j]])) / 
          length(comms[[i]])
      }
    }
    
    jaccard_mat <- jaccard_mat0
    diag(jaccard_mat) <- 0
    max_jacc <- max(jaccard_mat)
    deleted_comms <- integer(0)
    
    while (max_jacc > tau) {
      
      inds <- which(jaccard_mat == max_jacc, arr.ind = TRUE)[1, ]
      
      # keep comm with larger score
      delete_comm <- inds[which.min(c(scores[inds[1]], scores[inds[2]]))]
      jaccard_mat[delete_comm, ] <- 0
      jaccard_mat[, delete_comm] <- 0
      deleted_comms <- c(deleted_comms, delete_comm)
      max_jacc <- max(jaccard_mat)
      
    }
    
    kept_comms <- setdiff(1:K, deleted_comms)
    
    return(list("final_comms" = comms[kept_comms],
                "kept_comms" = kept_comms))
    
  }
  
  list2df <- function(listoflists) {
    #Usage: list2df(results$extraction_res)
    field_names <- unique(unlist(lapply(listoflists, names)))
    N <- length(listoflists)
    
    transpose_list <- rep(list(rep(list(), N)), length(field_names))
    names(transpose_list) <- field_names
    
    for(i in 1:N){
      for(n in field_names){
        transpose_list[[n]][i] <- listoflists[[i]][[n]] 
      }
    }
    
    do.call(data.frame, c(lapply(transpose_list, I), list(stringsAsFactors=FALSE)))
  }
  #-------------------------------------------------------------------------------
  # pvalue function
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
  #-------------------------------------------------------------------------------
  # BH procedure function
  bh_rejectR <- function (pvals, alpha, conserv = TRUE) {
    
    m <- length(pvals)
    
    if (!conserv) {
      pvals_adj <- m * pvals / rank(pvals, ties.method = "first")
    } else {
      mults <- sum(1 / c(1:m))
      pvals_adj <- mults * m * pvals / rank(pvals, ties.method = "first")
    }
    
    if (sum(pvals_adj <= alpha) > 0) {
      thres <- max(pvals[pvals_adj <= alpha])
      return(which(pvals <= thres))
    } else {
      return(integer(0))
    }
    
  }
  #-------------------------------------------------------------------------------
  # Initialize function
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
  #-------------------------------------------------------------------------------
  # Extraction function
  extract <- function (indx, interact = FALSE, print_output = verbose) {
    
    # Getting time of completion
    current_time <- proc.time()[3] - start_second
    
    # Doing checks
    if (interact && (current_time > time_limit || Dud_count > Dud_tol || OL_count > OL_tol)) {
      stop_extracting <<- TRUE
    }
    
    # If you want to interact with trackers, get the current tracking info
    if (interact) {
      
      # Checking for stops
      if (stop_extracting) return(list(report = "stop_extracting"))
      if (indx %in% clustered && !exhaustive) 
        return(list(report = "indx_clustered"))
      
    }
    
    # Print some output
    if (print_output && interact) {
      cat("\n#-----------------------------------------\n\n")
      cat("extraction", which(extractord == indx), "of", length(extractord), "\n\n")
      cat("indx", indx, "has not been clustered\n")
      cat("OL_count = ", OL_count, "\n")
      cat("Dud_count = ", Dud_count, "\n")
      cat(paste0(sum(clustered <= dx), " X vertices in communities.\n"))
      cat(paste0(sum(clustered > dx), " Y vertices in communities.\n"))
    }
    
    # Get general B01 & B02 (regardless of indx)
    B01 <- initialize1(indx, Cpp = Cpp)
    if (length(B01) > 1) {
      
      # Half-update
      pvals2 <- bmd_obj$pvals(B01, pval_parallel)
      if (Cpp) {
        B02 <- bh_rejectC(pvals2, alpha, conserv = TRUE)
      } else {
        B02 <- bh_rejectR(pvals2, alpha, conserv = TRUE)
      }
      
    } else {
      B02 <- integer(0)
    }
    
    # Check for dud
    if (length(B01) * length(B02) <= 1) {
      if (interact) {
        Dud_count <<- Dud_count + 1
      }
      if (updateOutput && verbose) {
        cat("--initialize was a dud\n")
      }
      return(list("indx" = indx, "report" = "dud"))
    }
    
    # Assign B0x and B0y according to indx
    if (indx > dx) {
      B02 <- Yindx[B02]
      B0x <- B01; B0y <- B02
    } else {
      B0x <- B02; B0y <- B01
    }
    
    if (updateOutput && verbose) {
      cat(paste0("--obtained initial set of size (", length(B0x), ", ", length(B0y), ")\n"))
    }
    
    # Initializing extraction loop
    B_oldx <- B0x; B_oldy <- B0y
    B_old <- c(B_oldx, B_oldy)
    initial_set <- B_old
    B_new <- c(Xindx, Yindx)
    chain <- list(B_old)
    consec_jaccards <- mean_jaccards <- NULL
    consec_sizes <- list(c(length(B_oldx), length(B_oldy)))
    found_cycle <- found_break <- cycledSets <- NULL
    did_it_cycle <- FALSE
    itCount <- 0
    
    # Extraction loop
    repeat {
      
      itCount <- itCount + 1
      cat(paste0("--calculating p-values for update ", itCount, "...\n"))
      
      if (updateMethod == 1) {
        
        pvalsXY <- c(bmd_obj$pvals(B_oldy, pval_parallel), bmd_obj$pvals(B_oldx, pval_parallel))
        B_new <- if (Cpp) {
          bh_rejectC(pvalsXY, alpha, conserv = TRUE)
        } else {
          bh_rejectR(pvalsXY, alpha)
        }
        B_newx <- B_new[B_new <= dx]; B_newy <- B_new[B_new > dx]
        
      }
      
      if (updateMethod == 2) {
        
        if (indx > dx) {
          
          # Update X nodes first
          pvalsx <- bmd_obj$pvals(B_oldy, pval_parallel)
          B_newx <- if(Cpp) {
            bh_rejectC(pvalsx, alpha, conserv = TRUE)
          } else {
            bh_rejectR(pvalsx, alpha, conserv = TRUE)
          }
          if (length(B_newx) <= 1) {
            B_newy <- integer(0)
            break
          }
          pvalsy <- bmd_obj$pvals(B_newx, pval_parallel)
          B_newy <- if (Cpp) {
            Yindx[bh_rejectC(pvalsy, alpha, conserv = TRUE)]
          } else {
            Yindx[bh_rejectR(pvalsy, alpha, conserv = TRUE)]
          }
          
        } else {
          
          # Update Y nodes first
          pvalsy <- bmd_obj$pvals(B_oldx, pval_parallel)
          B_newy <- if (Cpp) {
            Yindx[bh_rejectC(pvalsy, alpha, conserv = TRUE)]
          } else {
            Yindx[bh_rejectR(pvalsy, alpha, conserv = TRUE)]
          }
          if (length(B_newy) <= 1) {
            B_newx <- integer(0)
            break
          }
          pvalsx <- bmd_obj$pvals(B_newy, pval_parallel)
          B_newx <- if (Cpp) {
            bh_rejectC(pvalsx, alpha, conserv = TRUE)
          } else {
            bh_rejectR(pvalsx, alpha, conserv = TRUE)
          }
          
        }
        
        B_new <- c(B_newx, B_newy)
        
      }
      
      
      if (length(B_newy) * length(B_newx) == 0) {
        B_newy <- B_newx <- integer(0)
        if (updateOutput && verbose)
          cat("----collapse on at least one side\n")
        break
      }
      
      consec_jaccard <- jaccard(B_new, B_old)
      jaccards <- unlist(lapply(chain, function (B) jaccard(B, B_new)))
      found_cycle <- c(found_cycle, FALSE)
      found_break <- c(found_break, FALSE)
      consec_jaccards <- c(consec_jaccards, consec_jaccard)
      mean_jaccards <- c(mean_jaccards, mean(jaccards))
      consec_sizes <- c(consec_sizes, list(c(length(B_newx), length(B_newy))))
      
      # Checking for cycles (4.4.1 in CCME paper)
      if (jaccard(B_new, B_old) > 0) { # Otherwise loop will end naturally
        
        if (sum(jaccards == 0) > 0) { # Cycle has been found
          
          found_cycle[itCount] <- TRUE
          did_it_cycle <- TRUE
          
          if (updateOutput)
            cat("----cycle found")
          
          Start <- max(which(jaccards == 0))
          cycle_chain <- chain[Start:length(chain)]
          
          # Checking for cycle break (4.4.1a)
          seq_pair_jaccards <- rep(0, length(cycle_chain) - 1)
          for (j in seq_along(seq_pair_jaccards)) {
            seq_pair_jaccards[j] <- jaccard(chain[[Start + j - 1]], chain[[Start + j]])
          }
          if (sum(seq_pair_jaccards > 0.5) > 0) {# then break needed
            cat("------break found\n")
            found_break[itCount] <- TRUE
            B_new <- integer(0)
            break
          }
          
          # Create conglomerate set (and check, 4.4.1b)
          B_J <- unique(unlist(cycle_chain))
          B_J_check <- unlist(lapply(chain, function (B) jaccard(B_J, B)))
          if (sum(B_J_check == 0) > 0) {
            cat("------old cycle\n")
            break
          } else {
            cat("------new cycle\n")
            B_oldx <- B_J[B_J <= dx]
            B_oldy <- B_J[B_J > dx]
          }
          
        } else {
          # From checking jaccards to cycle_chain; if not, then can set B_oldx
          # and B_oldy to the update and restart.
          B_oldx <- B_newx
          B_oldy <- B_newy
        }
        
      } else { # From checking B_new to B_old; if not, B_new = B_old and:
        break
      }
      
      if (updateOutput && verbose) {
        cat(paste0("----updated to set of size (", length(B_newx), ", ", length(B_newy), ")\n"))
        cat(paste0("----jaccard to previous = ", round(consec_jaccard, 3), "\n"))
      }
      
      B_old <- c(B_oldx, B_oldy)
      chain <- c(chain, list(B_new))
      
    }
    
    # Storing B_new and collecting update info
    update_info <- list("mean_jaccards" = mean_jaccards, 
                        "consec_jaccards" = consec_jaccards,
                        "consec_sizes" = consec_sizes,
                        "found_cycle" = found_cycle,
                        "found_break" = found_break)
    if (length(B_newx) * length(B_newy) * length(B_new) == 0) {
      if (interact) {
        Dud_count <<- Dud_count + 1
      }
      return(list("indx" = indx,
                  "StableComm" = integer(0),
                  "update_info" = update_info,
                  "initial_set" = initial_set,
                  "itCount" = itCount, 
                  "did_it_cycle" = did_it_cycle,
                  "start_time" = current_time,
                  "report" = "break_or_collapse"))
    } else {
      clustered <<- union(clustered, B_new)
      comms <<- c(comms, list(B_new))
    }
    
    
    # Checking overlap with previous sets
    if (interact && length(comms) > 1) {
      OL_check <- unlist(lapply(comms, function (C) jaccard(B_new, C)))
      if (sum(OL_check < 1 - OL_thres) > 0) {
        OL_count <<- OL_count + 1
      }
    }
    
    return(list("indx" = indx,
                "StableComm" = B_new,
                "update_info" = update_info,
                "initial_set" = initial_set,
                "itCount" = itCount, "did_it_cycle" = did_it_cycle,
                "current_time" = current_time,
                "report" = "complete_extraction"))
    
  }
  
  #-------------------------------------------------------------------------------
  # Extractions
  
  if (generalOutput)
    cat("Beginning method.\n\n")
  
  # Getting node orders.
  Ysum <- Y %*% rep(1,dy) / dy
  Xsum <- X %*% rep(1,dx) / dx
  cor_X_to_Ysums <- abs(as.vector(t(Ysum) %*% X))
  cor_Y_to_Xsums <- abs(as.vector(t(Xsum) %*% Y))

  if (exhaustive) {
    extractord <- sample(c(Xindx, Yindx))
  } else {
    extractord <- c(Xindx, Yindx)[order(c(cor_X_to_Ysums, cor_Y_to_Xsums),
                                        decreasing = TRUE)]
  }
  
  if (!is.null(start_nodes))
    extractord <- extractord[extractord %in% start_nodes]
  
  # Initializing control variables
  stop_extracting <- FALSE
  OL_count <- 0
  Dud_count <- 0
  clustered <- integer(0)
  comms <- NULL
  
  # Extracting
  if (parallel) {
    ticp <- proc.time()[3]
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    clusterEvalQ(cl, library(bmdCpp))
    registerDoParallel(cl)
    extract_res <- foreach(i = extractord) %dopar% {
      extract(i, print_output = FALSE, interact = TRUE)
    }
    stopCluster(cl)
    tocp <- proc.time()[3]
  } else {
    tic <- proc.time()[3]
    extract_res <- lapply(extractord, extract, interact = TRUE)
    stopped <- sapply(extract_res, 
                      function (obj) obj$report == "stop_extracting")
    stopped_at <- min(which(stopped))
    extract_res <- extract_res[order(extractord)]
    toc <- proc.time()[3]
  }
  
  #-----------------------------------------------------------------------------
  # Clean-up and return -------------------------------------------------------
  
  if (generalOutput)
    cat("Cleaning up.\n")
  
  # Getting final sets and making report
  final.sets <- lapply(extract_res, function (R) R$StableComm)
  endtime <- proc.time()[3]
  report <- list(OL_count = OL_count, 
                 Dud_count = Dud_count, 
                 Total_time = endtime - start_second,
                 Extract_ord = extractord,
                 Stopped_at = stopped_at)
  names(report) <- c("OL_count", "Dud_count", "timer")
  
  # Removing blanks and trivial sets
  nonNullIndxs <- which(unlist(lapply(final.sets, length)) > 0)
  if (length(nonNullIndxs) == 0) {
    returnList <- list("communities" = list("X_sets" = NULL,
                                            "Y_sets" = NULL),
                       "background" = 1:(dx + dy),
                       "extract_res" = extract_res,
                       "finalIndxs" = integer(0),
                       "final.sets" = final.sets,
                       "report" = report)
    return(returnList)
  }
  nonNullComms <- final.sets[nonNullIndxs]
  OLfilt <- filter_overlap(nonNullComms, tau = OL_thres, inv.length = inv.length)
  finalComms <- OLfilt$final_comms
  finalIndxs <- nonNullIndxs[OLfilt$kept_comms]
  
  returnList <- list("extract_res" = extract_res,
                     "finalIndxs" = finalIndxs,
                     "report" = report)
  
  cat("#-------------------\n\n\n")
  
  return(returnList)
  
}
