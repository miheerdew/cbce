#' Correlation Bi-community Extraction method
#'
#' Given two groups of variables, find correlation bi-communities between them. For such a community, the nodes from the first group are are higly correlated to the community-nodes from the second group, and vice versa.
#'
#'  \code{cbce} applies an update function (mapping subsets of variables to subsets of variables) iteratively until a fixed point is found. These fixed points are reported as communities. The update function uses multiple-testing procedure to find variables correlated to the current set of variables.
#'
#' The update starts from a single node (starting with the initialization step) and is repeated till either a fixed point is found or some set repeats. Each such run is called an extraction. Since the extraction only starts from singleton node, there are \code{ncol(X)+ncol(Y)} possible extractions.
#'
#' @param X,Y Numeric Matices. Represents the two groups of variables.
#' @param alpha \eqn{\in (0,1)}. Controls the type1 error per update. This is the type1 error to use for multiple_testing procedure
#' @param alpha.init \eqn{\in (0,1)} Controls the type1 error for the initialization step. This should perhaps be more lax than alpha.
#' @param start_nodes The initial set of nodes to start with. If Null start from all the nodes.
#' @param calc_full_cor Calculate \code{c(ncol(X),ncol(Y))} dimensional correlation matrix. This makes the computation faster but requires more memory.
#' @param interaction This is a function that will be called between extracts to allow interaction with the program. For instance one cas pass interaction_gui (EXPERIMENTAL!) or interaction_none.  
#' @param multiple_testing_method Method to use for multiple testing. Should be one of c('BHY', 'BH', 'sqrt_BH', 'sqrt_BHY', 'square_BH', 'Bonferroni')
#' @param max_iterations The maximum number of iterations per extraction. If a fixed point is not found by this step, the extraciton is terminated.
#' @param diagnostic This is a function for probing the internal state of the method. It will be called at "Events" and can look into what the method is doing. Pass either diagnostics2, diagnostics_none or a custom function.
#' 
#' @return The return value is a list with details of the extraction and list of indices representing the communities. See example below (finding communities in noise). Note that the variables from the X and Y set are denoted using a single numbering. Hence the nodes in X are denoted by \code{1:dx} and the nodes in Y are denoted by the numbers following dx (hence \code{dx+1:dy})
#' 
#' @examples 
#' \dontrun{
#' n <- 100
#' dx <- 50
#' dy <- 70
#'
#' X <- matrix(rnorm(n*dx), ncol=dx)
#' Y <- matrix(rnorm(n*dy), ncol=dy)
#' res <- cbce2(X, Y)
#' df <- res$filtered_result.df
#' # The filtered bimodules:
#' bms <- rlist::list.map(res$extract_res[df$index], bimod)
#'}
#' @export
cbce2 <- function(X, Y, 
                  alpha = 0.05, 
                  alpha.init = alpha,
                  calc_full_cor=TRUE,
                  multiple_testing_method = 'BHY',
                  start_nodes=NULL,
                  max_iterations = 20,
                  interaction=interaction_none,
                  diagnostic=diagnostics2) {
  
  #--------------------------------
  # Defining useful functions.

  # Initialize the extraction.
  # Use the backend to do most of the work, but correct for global indices
  # @param indx The (global) index of the node (variable) to initialize from.
  # @return list(x=integer-vector, y=integer-vector): The initialized x, y sets.
  initialize <- function(indx) {
    B01 <- init(bk, indx, alpha.init, multiple_testing_method)
    if(length(B01) <= 1) {
      # If B01 = 1 declare it as dud.
      # This is because lot of B01 = 1 end up as duds
      # and we don't want to do the update step for them.
      return(list(x=integer(0), y=integer(0)))
    }
    if (indx <= dx) {
      #indx on the X side, so only need to correct the init-step.
      B01 <- B01 + dx
      B02 <- bh_reject(pvals(bk, B01), alpha.init, multiple_testing_method)
      return(list(x = B02, y = B01))
    } else {
      #indx on the Y side, so only need to correct the half update following the init step.
      B02 <- bh_reject(pvals(bk, B01), alpha.init, multiple_testing_method) + dx
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
  
  update <- function(B0, env=parent.frame(), first.update.X=TRUE) {
    # Do the update starting from B.
    # @param B0 list(x, y) : x is a subset of X nodes, y is a subset of Y nodes (using the global index).
    # @return list(x, y) : The X and Y subsets corresponding to the updated set (again using global numbering)
    
      B1 <- list()
      
      if(first.update.X) {
        B1$x <- bh_reject(pvals(bk, B0$y), alpha, multiple_testing_method)
        B1$y <- bh_reject(pvals(bk, B1$x), alpha, multiple_testing_method) + dx
      } else {
        B1$y <- bh_reject(pvals(bk, B0$x), alpha, multiple_testing_method) + dx
        B1$x <- bh_reject(pvals(bk, B1$y), alpha, multiple_testing_method)
      }
      
      #env$px <- px
      #env$py <- py
      
      #diagnostic("Update:Pvalues", env)
      #Note that the positions returned by bh_reject are already the global
      #numbers.
      return(B1)
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
      B1 <- update(B0, f, first.update.X = (indx > dx))
      
      if (length(B1$y) * length(B1$x) == 0) {
        stop <- collapsed <- TRUE
        diagnostic("Extract:Collapsed", f)
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
      } else if (did_it_cycle) {
        # Matches an old set.
        cycle_count = cycle_count + 1
        diagnostic("Extract:FoundCycle", f)
        if(cycle_count > 2) {
          # Cycled thrice already. Stop extraction.
          # This is because most iterations that cycle > 2 times
          # do not end up reaching a fixed point. And we don't have
          # infinite resources to keep dealing with cycles.
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
      log.pval <- summary_pval(X[, B0$x, drop=FALSE], Y[, B0$y, drop=FALSE])
    } else {
      stableComm <- list(x=integer(0), y=integer(0))
      log.pval <- NA
    }
    
    return(c(list("indx" = indx,
                  "bimod" = stableComm,
                  "itCount" = itCount,
                  "fixed_point" = success,
                  "collapsed" = collapsed,
                  "cycle_count" = cycle_count,
                  "overflowed" = !stop,
                  "log.pvalue" = log.pval), 
                  diagnostic_info))
  }
  
  #-------------------------------------------------------------------------------
  # Extractions
  
  
  # -------------- Global variables -------------------- 
  bk <- backend.perm(X, Y, calc_full_cor)
  
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
  list(extract_res=extract_res,
       filtered_res.df=filter_and_summarize(extract_res, show.progress = FALSE))
}