#' Correlation Bi-community Extraction method
#'
#' Consider two sets of high-dimensional measurements on the same 
#' set of samples. CBCE (Correlation Bi-Community Extraction method) 
#' finds sets of variables from the first measurement and sets of 
#' variables from the second measurement which are correlated to each 
#' other.
#'
#' \code{cbce} applies an update function (mapping subsets of
#'  variables to subsets of variables) iteratively until a fixed point 
#'  is found. These fixed points are reported as communities. 
#' The update starts from a single variable (the initialization step) 
#' and is repeated till either a fixed point is found or some set 
#' repeats. Each such run is called an extraction. Since the extraction
#'  only starts from singleton node, there are \code{ncol(X)+ncol(Y)} 
#'  possible extractions.
#'
#' @param X,Y Numeric Matices. Represents the two groups of variables. 
#' Rows represent samples and columns represent variables.
#' @param cov The covariates to account for; This should be a matrix 
#' with the same number of rows as X and Y. Each column represents
#' a covariate whose effect needs to be removed. If this is null,
#' no covariate will be removed.
#' @param alpha \eqn{\in (0,1)}. Controls the type1 error for the 
#' update (for the multiple testing procedure). 
#' @param alpha.init \eqn{\in (0,1)} Controls the type1 error 
#' for the initialization step. This could be more liberal 
#' (i.e greater than) than the alpha for the update step.
#' @param start_frac \eqn{\in (0,1)} The random proportion of nodes to 
#' start extractions from. This is used to randomly sample 
#' \code{start_nodes}. If \code{start_node} is provided this parameter 
#' is ignored.
#' @param start_nodes list The initial set of variables to start with.
#' If this is provided, \code{start_frac} will be ignored. 
#' If Null, extractions are run starting from each varable from X and Y.
#' Otherwise \code{start_node$x} gives the X variables to start from
#' and \code{start_nodes$y} gives the Y variables to start from.
#' @param cache.size integer The amount of memory to dedicate for 
#' caching correlations. This will speed things up. 
#' Defaults to the average memory required by X and Y matrices
#' @param max_iterations integer The maximum number of iterations per 
#' extraction. If a fixed point is not found by this step, 
#' the extraciton is terminated. This limit is set so that the
#' program terminates.
#' @param size_threshold The maximum size of bimodule we want to search for.
#'                    The search will be terminated when sets grow beyond this size.
#'                    The size of a bimodule is defined as the geometric mean of
#'                    its X and Y sizes.
#' @param interaction (internal) This is a function that will be called 
#' between extractions to allow interaction with the program. 
#' For instance one cas pass the function \code{\link{interaction_gui}}
#' (EXPERIMENTAL) or \code{\link{interaction_cli}}.  
#' @param diagnostic (internal) This is a internal function for 
#' probing the internal state of the method. It will be 
#' called at special hooks and can look into what the method is doing. 
#' Pass either \code{\link{diagnostics}}, 
#' \code{\link{diagnostics_none}}.
#' @param heuristic_search Use a fast, but incomplete, version of 
#' heuristic search that doesn't start from nodes inside bimodules 
#' already found.
#' @return The return value is a list with the results and 
#' meta-data about the extraction. The most useful field is
#' \code{comms} - this is a list of all the Correlation Bi-communities 
#' that was detected after filtering, while \code{comms.fil} 
#' consist of all the communities that were found after 
#' filtering similar communities.
#' 
#' @examples 
#' library(cbce)
#' #Sample size
#' n <- 40
#' #Dimension of measurement 1
#' dx <- 20
#' #Dimension of measurement 2
#' dy <- 50
#' #Correlation strength
#' rho <- 0.5
#' set.seed(1245)
#' # Assume first measurement is gaussian
#' X <- matrix(rnorm(dx*n), nrow=n, ncol=dx)
#' # Measurements 3:6 in set 2 are correlated to 4:7 in set 1
#' Y <- matrix(rnorm(dy*n), nrow=n, ncol=dy)
#' Y[, 3:6] <- sqrt(1-rho)*Y[, 3:6] + sqrt(rho)*rowSums(X[, 4:5])
#' res <- cbce(X, Y)
#' #Recovers the indices 4:5 for X and 3:6 for Y
#' #If the strength of the correlation was higher
#' #all the indices could be recovered.
#' res$comms
#' @export
cbce <- function(X, Y, 
                  alpha = 0.05, 
                  alpha.init = alpha,
                  cov = NULL,
                  cache.size = (utils::object.size(X) + 
                                utils::object.size(Y))/2,
                  start_frac = 1,
                  start_nodes = list(x=sample(1:ncol(X), 
                                            ceiling(ncol(X)*start_frac)), 
                                   y=sample(1:ncol(Y),
                                            ceiling(ncol(Y)*start_frac))),
                  max_iterations = 20,
                  size_threshold = 0.5*exp(log(ncol(X))/2 + log(ncol(Y))/2),
                  interaction=interaction_none,
                  heuristic_search=FALSE,
                  diagnostic=diagnostics) {
  
  #--------------------------------
  # Defining useful functions.

  # Initialize the extraction.
  # Use the backend to do most of the work, but correct for global 
  # indices
  # @param indx The (global) index of the node (variable) to 
  # initialize from.
  # @return list(x=integer-vector, y=integer-vector): The 
  # initialized x, y sets.
  initialize <- function(indx) {
    B01 <- rejectPvals(bk, indx, alpha.init)

    if (indx <= dx) {
      #indx on the X side, so only need to correct the init-step.
      B01 <- B01 + dx
      B02 <- rejectPvals(bk, B01, alpha.init)
      return(list(x = B02, y = B01))
    } else {
      #indx on the Y side, so only need to correct the half update 
      #following the init step.
      B02 <- rejectPvals(bk, B01, alpha.init) + dx
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
  # Note B$y is already uses global numbering, so no information is 
  # lost.
  merge <- function(B) {
    c(B$x, B$y)
  }
  
  update <- function(B0, env=parent.frame(), first.update.X=TRUE) {
    # Do the update starting from B.
    # @param B0 list(x, y) : x is a subset of X nodes, y is a subset 
    # of Y nodes (using the global index).
    # @return list(x, y) : The X and Y subsets corresponding to the 
    #updated set (again using global numbering)
    
      B1 <- list()
      
      if(first.update.X) {
        B1$x <- rejectPvals(bk, B0$y, alpha)
        B1$y <- rejectPvals(bk, B1$x, alpha) + dx
      } else {
        B1$y <- rejectPvals(bk, B0$x, alpha) + dx
        B1$x <- rejectPvals(bk, B1$y, alpha)
      }
      
      #env$px <- px
      #env$py <- py
      
      #diagnostic("Update:Pvalues", env)
      #Note that the positions returned by bh_reject are already the 
      #global
      #numbers.
      return(B1)
  }
  
  hash <- function(B) {
    digest::digest(c(sort(B$x), sort(B$y)))
  }
  
  extract <- function(indx) {
    # Start the extraction from indx
    #
    # First do the initialization at indx, and then repeateadly apply 
    #update till a fixed point is found, or one of the sets repeates 
    #(forming a non-trivial cycle).
    #@return The return value is the extract_res field of the final
    #method results.
    
    #Is indx an X or Y node. This is important for the two-step update
    #to tell \code{update()} which side to start initialization from.
    

    cycle_count <- 0
    collapsed <- FALSE
    success <- FALSE
    itCount <- 0
    size_exceeded <- FALSE
    stop <- FALSE

    f <- new.env()
    diagnostic("Extract:Setup", f)
    
    B0 <- initialize(indx)
    
    # Check for dud
    if (length(B0$x)*length(B0$y) <= 0) {
      collapsed <- stop <- TRUE
      #diagnostic('Extract:Dud', f)
    }
    
    # Initializing extraction loop
    
    #The hashes of all sets visited during the extraction
    chain <- rep(NULL, max_iterations) 
    
    B1 <- split(integer(0)) #Initial value for B1
    
    diagnostic("Extract:LoopBegins", f)
    # Extraction loop
    while(itCount < max_iterations && !stop) {
      itCount <- itCount + 1
      B1 <- update(B0, f, first.update.X = (indx > dx))
      
      bm.size <- sqrt(length(B1$y)*length(B1$x))
      
      if (bm.size < 1) {
        stop <- collapsed <- TRUE
        diagnostic("Extract:Collapsed", f)
        break
      }
      
      if(bm.size > size_threshold) {
        size_exceeded <- stop <- TRUE
        diagnostic("Extract:SizeExceeded", f)
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
      log.pval <- summary_pval(X[, B0$x, drop=FALSE], 
                               Y[, B0$y, drop=FALSE],
                               n.eff=n.eff)
    } else {
      stableComm <- list(x=integer(0), y=integer(0))
      log.pval <- NA
    }
    
    return(c(list(
                  #Convert index to local number 
                  "indx" = if(indx > dx) c(y=indx-dx) else c(x=indx),
                  "bimod" = stableComm,
                  "itCount" = itCount,
                  "fixed_point" = success,
                  "collapsed" = collapsed,
                  "cycle_count" = cycle_count,
                  "overflowed" = !stop,
                  "size_exceeded" = size_exceeded,
                  "log.pvalue" = log.pval), 
                  diagnostic_info))
  }
  
  #-----------------------------------------------------
  # Extractions
  
  if(!is.null(cov)) {
    # ------------- Covariate correction -----------------
    # center the covariance matrix
    cov <- scale(cov, scale=FALSE)
    q <- qr(cov)
    Q <- qr.Q(q)
    
    #Residualize the X and Y matrices
    X <- X - Q %*% crossprod(Q,X)
    Y <- Y - Q %*% crossprod(Q,Y)
    #Reduce the effective dimension effective dimension
    n.eff <- nrow(X) - q$rank
  } else {
    n.eff <- nrow(X)
  }
  # -------------- Global variables --------------------
  bk <- backend.perm(X, Y, cache.size, n.eff=n.eff)
  
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
  
  if (!is.null(start_nodes)) {
    start_nodes_global <- c(start_nodes$x, start_nodes$y + dx)
    if(is.null(start_nodes_global)) {
      warning(paste0("start_nodes must be a list of size 2",
                  "and must have names x and y. Ignoring it."))
    } else {
      extractord <- extractord[extractord %in% start_nodes_global]
    }
  }
  # --------- Extraction loop starts ----------------
  extract_res <- rlist::list.map(1:length(extractord), NULL)
  res <- NULL
  
  #Create a new env to pass to interaction() 
  e <- new.env()
  
  # All the nodes in the bimodules found so far.
  # Helps with heuristic search.
  bimod.nodes <- hash::hash()
  
  interaction("Main:Setup", e)
  
  for(i in seq_along(extractord)) {
    if(heuristic_search && 
        hash::has.key(hash::make.keys(i), bimod.nodes)
       ) {
      next
    }
    
    res <- extract(extractord[i])
    
    if(res$fixed_point) {
      bimod.nodes[hash::make.keys(
        c(res$bimod$x, res$bimod$y + dx))] <- TRUE
    }
    
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
  
  extract_res <- rlist::list.filter(extract_res, !is.null(.))
    
  interaction("Main:Filtering", e)
  summary <- filter_and_summarize(extract_res, logpval.thresh=log(alpha)-(log(dx)+log(dy)))
  
  interaction("Main:End", e)
  
  list(extract_res=extract_res,
       summary=summary,
       comms.fil=rlist::list.map(extract_res[summary$df.fil$index], bimod),
       comms=rlist::list.map(extract_res[summary$df.unique$index],
                                 bimod)
  )
}
