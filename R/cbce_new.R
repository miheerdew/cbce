#' @export
cbce2 <- function(X, Y, 
                  alpha = 0.05, 
                  init_method = "conservative-BH",
                  calc_full_cor=TRUE,
                  multiple_testing_method = 'BHY',
                  start_nodes=NULL,
                  max_iterations = 20,
                  break_thresh = 0.6, #If the disjointness > break_tresh, then stop iteration.
                  interaction=interaction_none,
                  diagnostic=diagnostics2) {
  
  #--------------------------------
  # Defining useful functions.

  # Initialize the extraction.
  # Use the backend to do most of the work, but correct for global indices
  # @param indx The (global) index of the node (variable) to initialize from.
  # @return list(x=integer-vector, y=integer-vector): The initialized x, y sets.
  initialize <- function(indx) {
    B01 <- init(bk, indx, alpha, init_method)
    if(length(B01) <= 1) {
      # If B01 = 1 declare it as dud.
      # This is because lot of B01 = 1 end up as duds
      # and we don't want to do the update step for them.
      return(list(x=integer(0), y=integer(0)))
    }
    if (indx <= dx) {
      #indx on the X side, so only need to correct the init-step.
      B01 <- B01 + dx
      B02 <- bh_reject(pvals(bk, B01), alpha, multiple_testing_method)
      return(list(x = B02, y = B01))
    } else {
      #indx on the Y side, so only need to correct the half update following the init step.
      B02 <- bh_reject(pvals(bk, B01), alpha, multiple_testing_method) + dx
      return(list(x = B01, y = B02))
    }
  }
  
  #Split a single set into X, Y each with the global numbering.
  split <- function(set){
    set_x <- Filter(function(x) x <= dx, set)
    set_y <- Filter(function(x) x > dx, set)
    return(list(x=set_x, y=set_y))
  }
  
  # Store the X and Y sets of B in a single vector.
  # This is an inverse to split
  # Note B$y is already uses global numbering, so no information is lost.
  merge <- function(B) {
    c(B$x, B$y)
  }
  
  update <- function(B0, env=parent.frame()) {
    # Do the update starting from B.
    # @param B0 list(x, y) : x is a subset of X nodes, y is a subset of Y nodes (using the global index).
    # @return list(x, y) : The X and Y subsets corresponding to the updated set (again using global numbering)
    
      px <- pvals(bk, B0$y) # size dx vector
      py <- pvals(bk, B0$x) # size dy vector
      
      env$px <- px
      env$py <- py
      
      diagnostic("Update:Pvalues", env)
      #Note that the positions returned by bh_reject are already the global
      #numbers.
      return(split(bh_reject(c(px, py), alpha, multiple_testing_method)))
  }
  
  hash <- function(B) {
    digest::digest(c(sort(B$x), sort(B$y)))
  }
  
  extract <- function(indx) {
    # Start the extraction from indx
    #
    # First do the initialization at indx, and then repeateadly apply update till a fixed point is found, or one of the sets repeates (forming a non-trivial cycle).
    #
    #@return The return value is the extract_res field of the final method results.
    
    #Is indx an X or Y node. This is important for the two-step update to
    # tell \code{update()} which side to start initialization from.
    

    cycle_count <- 0
    collapsed <- FALSE
    success <- FALSE
    itCount <- 0
    stop <- FALSE

    f <- new.env()
    diagnostic("Extract:Setup", f)
    
    B0 <- initialize(indx)
    
    # Check for dud
    if (length(B0$x)*length(B0$y) <= 1) {
      collapsed <- stop <- TRUE
      #diagnostic('Extract:Dud', f)
    }
    
    # Initializing extraction loop
    
    #The hashes of all sets visited during the extraction
    chain <- rep(NULL, max_iterations) 
    
    B1 <- split(integer(0)) #Initial value for B1
    
    
    # Extraction loop
    while(itCount < max_iterations && !stop) {
      itCount <- itCount + 1
      B1 <- update(B0, f)
      
      if (length(B1$y) * length(B1$x) == 0) {
        stop <- collapsed <- TRUE
        #diagnostic("Extract:Collapsed", f)
        break
      }
      
      dist_to_prev <- jaccard_pairs(B1, B0)
      h <- hash(B1)
      did_it_cycle <- h %in% chain
      diagnostic("Extract:AfterUpdate", f)
      
      chain[itCount] <- h
      
      
      if (dist_to_prev == 0) {
        #End loop if fixed point found.
        stop <- success <- TRUE
      } else if (disjointness_pairs(B1, B0) > break_thresh) {
        #The sets are too disjoint. Start from the intersection.
        diagnostic("Extract:Disjoint", f)
        B0 <- intersect_pairs(B0, B1)
      } else if (did_it_cycle) {
        # Matches an old set.
        cycle_count = cycle_count + 1
        diagnostic("Extract:FoundCycle", f)
        if(cycle_count > 2) {
          #Cycled thrice already. Stop extraction.
          diagnostic("Extract:FoundBreak", f)
          stop <- TRUE
          B0 <- B1
        } else {
          # After a cycle, give an opportunity to start 
          # from the intersection.
          B0 <- intersect_pairs(B0, B1)
        }
      } else {
        B0 <- B1
      }
    }
    
    diagnostic_info <- diagnostic("Extract:End", f)

    if(!collapsed) {
      B0$y <- B0$y - dx
      stableComm <- B0 
    } else {
      stableComm <- list(x=integer(0), y=integer(0))
    }
    
    return(c(list("indx" = indx,
                  "bimod" = stableComm,
                  "itCount" = itCount,
                  "fixed_point" = success,
                  "collapsed" = collapsed,
                  "cycle_count" = cycle_count,
                  "overflowed" = !stop), 
                  diagnostic_info))
  }
  
  #-------------------------------------------------------------------------------
  # Extractions
  
  
  # -------------- Global variables -------------------- 
  bk <- backend.chisq(X, Y, calc_full_cor)
  
  dx <- ncol(X)
  dy <- ncol(Y)
  # ----------------------------------------------------
  
  indices <- c(1:dx, 1:dy + dx)
  
  # Getting node orders. We know what nodes to start with, 
  # but what is a good order to start in?
  

  #Use the sum of correlation hurestic to decide the order
  #bk$X, bk$Y are scaled versions of X and Y
  Ysum <- bk$Y %*% rep(1,dy) / dy
  Xsum <- bk$X %*% rep(1,dx) / dx
  cor_X_to_Ysums <- abs(as.vector(t(Ysum) %*% bk$X))
  cor_Y_to_Xsums <- abs(as.vector(t(Xsum) %*% bk$Y))
  extractord <- indices[order(c(cor_X_to_Ysums, cor_Y_to_Xsums),
                                decreasing = TRUE)]
  if (!is.null(start_nodes))
    extractord <- extractord[extractord %in% start_nodes]
  
  # --------- Extraction loop starts ----------------
  extract_res <- rlist::list.map(1:length(extractord), NULL)
  res <- NULL
  
  #Create a new env to pass to interaction() 
  e <- new.env()
  
  interaction("Main:Setup", e)
  
  for(i in seq_along(extractord)) {
    res <- extract(extractord[i])
    
    action <- interaction("Main:NextExtraction", e)
    
    if(!is.null(action) && action %in% c("stop", "browse")) {
      if(action == "browse") {
        browser()
      } else {
        break
      }
    } 
     
    extract_res[[i]] <- res
  }
  
  interaction("Main:End", e)
  list(extract_res=extract_res)
}