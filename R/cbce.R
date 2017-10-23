#source("backend.R")
#source("diagnostics.R")
#source("helper.R")

#' Correlation Bi-community Extraction method 
#' 
#' Given two groups of variables, find correlation bi-communities between them. \code{cbce} applies an update function (mapping subsets of variables to subsets of variables) iteratively until a fixed point is found. These fixed points are reported as communities. The update function uses multiple-testing procedure to find variables correlated to the current set of variables.
#' 
#' The update starts from a single node (starting the initialization step) and then is repeated till either a fixed point is found or some set repeats. Each such run is called an extraction. Since the extraction only starts from singleton node, there are \code{ncol(X)+ncol(Y)} possible extractions.
#' 
#' @param X,Y Numeric Matices. Represents the two groups of variables.
#' @param alpha Numeric (default 0.05). Controls the type1 error per update.
#' @param exhaustive Boolean (default FALSE). If exhaustive is FALSE, new extractions are not started from nodes within found communities. Otherwise, attempt is made to start from all communities.
#' @return The return value is a list with details of the extraction and list of indices representing the communities. See example below (finding communities in noise). Note that the variables from the X and Y set are denoted using a single numbering. Hence the nodes in X are denoted by 1:dx and the nodes in Y are denoted by the numbers following dx (hence dx+1:dy)
#' @example 
#' n <- 100
#' dx <- 200
#' dy <- 100
#' 
#' X <- matrix(rnorm(n*dx), ncol=dx)
#' Y <- matrix(rnorm(n*dy), ncol=dy)
#' res <- cbce(X, Y)
#' finalComms <- lapply(result$extract_res[res$finalIndxs], function(r) r$StableComm)
cbce <- function(
                X, Y,  #The X and Y matrices
                alpha = 0.05,  # The type1 error to use for BH procedure
                OL_thres = 0.9, # Thresholds used for 
                exhaustive = FALSE, #Initialize starting from all Gene/SNP nodes
                OL_tol = Inf, #If more than OL_tol found communities have OL_thresh overlap, stop method.
                Dud_tol = Inf, #If more than Dud_tol initializations end in a dud, stop method.
                time_limit = 18000, #Stop method after time_limit seconds.
                updateMethod = 1, #use the 1-step update vs 2-step update
                init_method = 1,
                inv.length = TRUE, #?use inv.length while disjointifying 
                start_nodes = NULL, #The initial set of nodes to start with
                parallel = FALSE, # Do the p-value computation in parallel
                calc_full_cor=FALSE, # Calculate ncol(X) \times ncol(Y) correlation matrix
                backend = "chisq", #One in "normal", "normal_two_sided", "chisq"
                diagnostics=diagnostics1
                ) {
  
  #-----------------------------------------------------------------------------
  # Setup 
  
  start_second <- proc.time()[3]
  cat("#-------------------\n")
  
  cat("Setting up backend\n")
  #Initialize the backend method specified in \code{backend}. This involves doing some precomputation. The variable bk is an S3 object, it holds the precomputations necessary to do the pvalue computations.
  bk <- getS3method("backend", backend)(X, Y, calc_full_cor, init_method)
   
  dx <- ncol(X) 
  dy <- ncol(Y)
  
  #Global numbering scheme for the X and Y nodes.
  #Remember how X indices were numbered from (1:dx) and Y from (dx+1 : dy)?
  Xindx <- 1:dx 
  Yindx <- (1:dy) + dx
  
  
  cat("Beginning method.\n\n")
  
  # Getting node orders. We know what nodes to start with, but what is a good order to start in.
  if (exhaustive) {
   # Random ordering
    extractord <- sample(c(Xindx, Yindx))
  } else {
    # Sort the nodes in decreasing order of thier correlation sums to the oppisite side.
    # This can be calculated efficiently. Note that bk$X, bk$Y are scaled verions of X and Y.
    Ysum <- bk$Y %*% rep(1,dy) / dy
    Xsum <- bk$X %*% rep(1,dx) / dx
    cor_X_to_Ysums <- abs(as.vector(t(Ysum) %*% bk$X))
    cor_Y_to_Xsums <- abs(as.vector(t(Xsum) %*% bk$Y))
    extractord <- c(Xindx, Yindx)[order(c(cor_X_to_Ysums, cor_Y_to_Xsums),
                                        decreasing = TRUE)]
  }
  
  if (!is.null(start_nodes))
    extractord <- extractord[extractord %in% start_nodes]
  
  # Initializing control variables
  stop_extracting <- FALSE
  # The Dud and Overlap counts
  OL_count <- 0
  Dud_count <- 0
  # The collection of nodes already clustered.
  clustered <- integer(0)
  comms <- NULL
  
  #-------------------------------------------------------------------------------
  # The core functions
  
  
  #' Initialize the extraction.
  #' Use the backend to do most of the work, but correct for global indices
  #' @param indx The (global) index of the node (variable) to initialize from.
  #' @return list(x=integer-vector, y=integer-vector): The initialized x, y sets. 
  initialize <- function(indx) {
    if (indx <= dx) {
      #indx on the X side, so only need to correct the init-step.
      B01 <- init(bk, indx, alpha) + dx
      if(length(B01) > 1) {
        B02 <- bh_reject(pvals(bk, B01), alpha)
      } else {
        B02 <- integer(0)
      }
      return(list(x = B02, y = B01))
    } else {
      #indx on the Y side, so only need to correct the half update following the init step.
      B01 <- init(bk, indx, alpha)
      if(length(B01) > 1) {
        B02 <- bh_reject(pvals(bk, B01), alpha) + dx
      } else {
        B02 <- integer(0)
      }
      return(list(x = B01, y = B02))
    }
  }

  #Split a single set into X, Y each with the global numbering.
  split_half <- function(set, dx){
    set_x <- Filter(function(x) x <= dx, set)
    set_y <- Filter(function(x) x > dx, set)
    return(list(x=set_x, y=set_y))
  } 

  update <- function(B0, startX=TRUE) {
  #' Do the update starting from B. Do either the two sided or one-sided update.
  #' @param B0 list(x, y) : x is a subset of X nodes, y is a subset of Y nodes (using the global index).
  #' @param startX If doing 2 step update start at X side; If FALSE,
  #'  start at Y side. Does not have any effect on 1 step update.
  #' @return list(x, y) : The X and Y subsets corresponding to the updated set (again using global numbering) 
    
    if (updateMethod == 2) {
      if (startX) { 
        B1y <- bh_reject(pvals(bk, B0$x), alpha) + dx
        B1x <- bh_reject(pvals(bk, B1y), alpha)
      } else {
        B1x <- bh_reject(pvals(bk, B0$y), alpha)
        B1y <- bh_reject(pvals(bk, B1x), alpha) + dx
      }
      return(c(B1x,B1y))
    } else {
      px <- pvals(bk, B0$y) # size dx vector
      py <- pvals(bk, B0$x) # size dy vector
      
      #Note that the positions returned by bh_reject are already the global 
      #numbers.
      return(split_half(bh_reject(c(px, py), alpha), dx))
    }
  }
  
  extract <- function(indx) {
    #' Start the extraction from indx
    #' 
    #' First do the initialization at indx, and then repeateadly apply update till a fixed point is found, or one of the sets repeates (forming a non-trivial cycle). 
    #' In case there is a cycle and if none of the consecutive sets in the cycle are sufficiently disjoint (jaccard >= 0.5), the take the union of all those sets and start the extraction from the conglomerate set. Othewise we say is break is found, and the extraction is terminated (without any result).
    #' 
    #'@return The return value is the extract_res field of the final method results.
    
    # Getting start time
    current_time <- proc.time()[3] - start_second
    
    # Doing checks
    if (current_time > time_limit || Dud_count > Dud_tol || OL_count > OL_tol) {
      stop_extracting <<- TRUE
    }
    
    if (stop_extracting) return(list(report = "stop_extracting"))
    
    # Checking for stops
    if (indx %in% clustered && !exhaustive) 
      return(list(report = "indx_clustered"))
    
    # Print some output
    cat("\n#-----------------------------------------\n\n")
    cat("extraction", which(extractord == indx), "of", length(extractord), "\n\n")
    cat("indx", indx, "has not been clustered\n")
    cat("OL_count = ", OL_count, "\n")
    cat("Dud_count = ", Dud_count, "\n")
    cat(paste0(sum(clustered <= dx), " X vertices in communities.\n"))
    cat(paste0(sum(clustered > dx), " Y vertices in communities.\n"))
    
    #Is indx an X or Y node. This is important for the two-step update to 
    # tell \code{update()} which side to start initialization from.
    startsAtX <- (indx <= dx) 
    
    B0 <- initialize(indx)
    
    # Check for dud
    if (length(B0$x) * length(B0$y) <= 1) {
      Dud_count <<- Dud_count + 1
      cat("--initialize was a dud\n")
      return(list("indx" = indx, "report" = "dud"))
    }

    cat(paste0("--obtained initial set of size (", length(B0$x), ", ", length(B0$y), ")\n"))
    
    # Initializing extraction loop
    
    B_old <- c(B0$x, B0$y)  #Store the starting set in a single vector
    
    chain <- list(B_old) #The collection of all sets visited during the extraction
    
    B_new <- c(Xindx, Yindx)
    did_it_cycle <- FALSE
    itCount <- 0
    
    diagnostics("ExtracionLoopBegins")
    
    # Extraction loop
    repeat {
      itCount <- itCount + 1
      cat(paste0("--Performing update ", itCount, "...\n"))
      
      B1 <- update(B0, startX = startsAtX)
      
      if (length(B1$y) * length(B1$x) == 0) {
        B1$y <- B1$x <- integer(0)
        diagnostics("Collapsed")
        cat("----collapse on at least one side\n")
        break
      }
      
      # Again, store B0, B1 in a single vector
      B_new <- unlist(B1, use.names = FALSE)
      B_old <- unlist(B0, use.names = FALSE)
      
      jaccards <- list.mapv(chain, jaccard(., B_new))

      diagnostics("AfterUpdate")
      
      # Checking for cycles (4.4.1 in CCME paper)
      if (jaccard(B_new, B_old) > 0) { # Otherwise loop will end naturally
        
        if (sum(jaccards == 0) > 0) { # Cycle has been found
          diagnostics("FoundCycle")

          did_it_cycle <- TRUE
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
            diagnostics("FoundBreak")
            B_new <- integer(0)
            break
          }
          
          # Create conglomerate set (and check, 4.4.1b)
          B_J <- unique(unlist(cycle_chain))
          B_J_check <- list.mapv(chain, jaccard(B_J, .))
          if (sum(B_J_check == 0) > 0) {
            cat("------old cycle\n")
            break
          } else {
            cat("------new cycle\n")
            B0 <- split_half(B_J, dx)
          }
          
        } else {
          # From checking jaccards to cycle_chain; if not, then can set B0
          # to the update and restart.
          B0 <- B1
        }
        
      } else {# From checking B_new to B_old; if not, B_new = B_old and:
        break
      }
      cat(paste0("----updated to set of size (", length(B1$x), ", ", length(B1$y), ")\n"))
      cat(paste0("----jaccard to previous = ", round(jaccards[length(jaccards)], 3), "\n"))
      chain <- c(chain, list(B_new))
    }
    
    diagnostic_info <- diagnostics("EndOfExtract")
    # Storing B_new and collecting update info
  
    if (length(B_new) == 0) {
     Dud_count <<- Dud_count + 1
     return(c(list("indx" = indx,
                 "StableComm" = integer(0),
                 "itCount" = itCount, 
                 "did_it_cycle" = did_it_cycle,
                 "report" = "break_or_collapse"), diagnostic_info))
    } else {
      clustered <<- union(clustered, B_new)
      comms <<- list.append(comms, B_new)
    }
    
    
    # Checking overlap with previous sets
    if (length(comms) > 1) {
      OL_check <- list.mapv(comms, jaccard(B_new, .))
      if (sum(OL_check < 1 - OL_thres) > 0) {
        OL_count <<- OL_count + 1
      }
    }
    
    return(c(list("indx" = indx,
                "StableComm" = B_new,
                "itCount" = itCount, "did_it_cycle" = did_it_cycle,
                "report" = "complete_extraction"), diagnostic_info))
    
  }
  
  #-------------------------------------------------------------------------------
  # Extractions
  
  
  # Extracting
  extract_res <- lapply(extractord, extract)
  stopped_at <-  list.first(extract_res, report == "stop_extracting")
  extract_res <- extract_res[order(extractord)]
  
  #-----------------------------------------------------------------------------
  # Clean-up and return -------------------------------------------------------
  
  cat("Cleaning up.\n")
  
  # Getting final sets and making report
  final.sets <- lapply(extract_res, function(R) R$StableComm)
  endtime <- proc.time()[3]
  report <- list(OL_count = OL_count, 
                 Dud_count = Dud_count, 
                 Total_time = endtime - start_second,
                 Extract_ord = extractord,
                 Stopped_at = stopped_at)
  
  # Removing blanks and trivial sets
  nonNullIndxs <- list.which(final.sets, length(.) > 0)
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
  # Keep a single representative for communities with overlap > OL_thres
  OLfilt <- filter_overlap(nonNullComms, tau = OL_thres, inv.length = inv.length)
  #finalComms <- OLfilt$final_comms
  finalIndxs <- nonNullIndxs[OLfilt$kept_comms]
  
  returnList <- list("extract_res" = extract_res,
                     "finalIndxs" = finalIndxs,
                     "report" = report)
  
  cat("#-------------------\n\n\n")
  
  return(returnList)
}
