#' Correlation Bi-community Extraction method
#'
#' Given two groups of variables, find correlation bi-communities between them. For such a community, the nodes from the first group are are higly correlated to the community-nodes from the second group, and vice versa.
#'
#'  \code{cbce} applies an update function (mapping subsets of variables to subsets of variables) iteratively until a fixed point is found. These fixed points are reported as communities. The update function uses multiple-testing procedure to find variables correlated to the current set of variables.
#'
#' The update starts from a single node (starting with the initialization step) and is repeated till either a fixed point is found or some set repeats. Each such run is called an extraction. Since the extraction only starts from singleton node, there are \code{ncol(X)+ncol(Y)} possible extractions.
#'
#' @param X,Y Numeric Matices. Represents the two groups of variables.
#' @param alpha \eqn{\in (0,1)}. Controls the type1 error per update. This is the type1 error to use for BH procedure
#' @param exhaustive Boolean. If exhaustive is FALSE, new extractions are not started from nodes within found communities. Otherwise, attempt is made to start from all nodes.
#' @param OL_thres \eqn{\in (0,1)}. Threshold  used to conclude significant overlap. The sets have significant overlap if jaccard similarity is >= OL_thres.
#' @param OL_tol If more than OL_tol found communities have OL_thres overlap, stop method.
#' @param Dud_tol If more than Dud_tol initializations end in a dud, stop method.
#' @param time_limit Stop method after time_limit seconds.
#' @param updateMethod Use the 1(-step) update vs 2(-step) update
#' @param init_method The initialization procedure to use. Must be one of "conservative-BH", "non-conservative-BH", "BH-0.5",  "BH-0.5-nc", "BH-0.9-nc", "no-multiple-testing".
#' @param inv.length Logical. Use inv.length as score while selecting the smallest community from communities with significant overlap.
#' @param start_nodes The initial set of nodes to start with. If Null start from all the nodes (may still exclude nodes within found communities if exhaustive = FALSE).
#' @param parallel Use parallel processing.
#' @param calc_full_cor Calculate \code{c(ncol(X),ncol(Y))} dimensional correlation matrix. This makes the computation faster but requires more memory.
#' @param mask_extracted Boolean. If true, then once a stable community is found, the nodes in it will never be used in future extractions.
#' @param backend The engine to use for p-value computation. Currently must be one of "normal", "normal_two_sided", "chisq", "chisq_fast", "indepChiSq", "indepNormal".
#' @param diagnostics This is a function that is called whenever internal events happen. It can then collect useful meta-data which is added to the final resutls.
#' @param rank_initial_sets Logical Start extraction from inital sets with the higher scores. 
#' @param init_quick_update Logical Use a quick half-update in the init step.
#' @param break_thresh Numeric A number between [0,1] that determines when a cycle will be broken. A value of 1 means that cycles will allways be broken and 0 that it will never be.
#' @return The return value is a list with details of the extraction and list of indices representing the communities. See example below (finding communities in noise). Note that the variables from the X and Y set are denoted using a single numbering. Hence the nodes in X are denoted by \code{1:dx} and the nodes in Y are denoted by the numbers following dx (hence \code{dx+1:dy})
#' @export
#' @examples
#' \dontrun{
#' n <- 100
#' dx <- 50
#' dy <- 70
#'
#' X <- matrix(rnorm(n*dx), ncol=dx)
#' Y <- matrix(rnorm(n*dy), ncol=dy)
#' res <- cbce(X, Y)
#' finalComms <- lapply(res$extract_res[res$finalIndxs], function(r) r$StableComm)
#'}
cbce <- function(X, Y,
                alpha = 0.05,
                OL_thres = 0.9,
                break_thresh = 0.5,
                exhaustive = FALSE,
                OL_tol = Inf,
                Dud_tol = Inf,
                time_limit = Inf,
                updateMethod = 1,
                init_method = "conservative-BH",
                inv.length = TRUE,
                start_nodes = NULL,
                parallel = FALSE,
                calc_full_cor=FALSE,
                backend = "chisq",
                mask_extracted = FALSE,
                rank_initial_sets = FALSE,
                init_quick_update = FALSE,
                diagnostics=diagnostics1
                ) {

  #-----------------------------------------------------------------------------
  # Setup

  if(mask_extracted && exhaustive) {
    stop("mask_extracted and exhaustive are incompatible.")
  }

  start_second <- proc.time()[3]
  cat("#-------------------\n")

  cat("Setting up backend\n")
  #Initialize the backend method specified by \code{backend}. This involves doing some precomputation. The variable bk is an S3 object which stores the precomputation.
  bk <- switch(backend,
               chisq = backend.chisq(X, Y, parallel, calc_full_cor),
               chisq_fast = backend.chisq(X, Y, parallel, calc_full_cor, fast_approx=TRUE),
               normal = backend.normal(X, Y, calc_full_cor),
               normal_two_sided = backend.normal_two_sided(X, Y, calc_full_cor),
               indepChiSq = backend.indepChiSq(X, Y, calc_full_cor),
               indepNormal = backend.indepNormal(X, Y, calc_full_cor),
               stop(paste("Unknown backend:", backend))
      )

  dx <- ncol(X)
  dy <- ncol(Y)

  #Global numbering scheme for the X and Y nodes.
  #Remember how X indices were numbered from (1:dx) and Y from (dx+1 : dy)?
  Xindx <- 1:dx
  Yindx <- (1:dy) + dx


  cat("Beginning method.\n\n")

  # Initializing control variables
  stop_extracting <- FALSE
  # The Dud and Overlap counts
  OL_count <- 0
  Dud_count <- 0
  # The collection of nodes already clustered.
  clustered <- integer(0)
  comms <- list()
  pvals_init <- if (init_quick_update) pvals_quick else pvals

  #-------------------------------------------------------------------------------
  # The core functions


  # Initialize the extraction.
  # Use the backend to do most of the work, but correct for global indices
  # @param indx The (global) index of the node (variable) to initialize from.
  # @return list(x=integer-vector, y=integer-vector): The initialized x, y sets.
  initialize <- function(indx, pval_func=pvals) {
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
      B02 <- bh_reject(pval_func(bk, B01), alpha)
      return(list(x = B02, y = B01))
    } else {
      #indx on the Y side, so only need to correct the half update following the init step.
      B02 <- bh_reject(pval_func(bk, B01), alpha) + dx
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

  jac <- function(B1, B2) {
    jaccard(merge(B1), merge(B2))
  }

  conglomerate <- function(chain) {
    #Construct the union of the sets in chain
    split(unique(unlist(chain)))
  }

  update <- function(B0, startX=TRUE) {
  # Do the update starting from B. Do either the two sided or one-sided update.
  # @param B0 list(x, y) : x is a subset of X nodes, y is a subset of Y nodes (using the global index).
  # @param startX If doing 2 step update start at X side; If FALSE,
  #  start at Y side. Does not have any effect on 1 step update.
  # @return list(x, y) : The X and Y subsets corresponding to the updated set (again using global numbering)

    if (updateMethod == 2) {
      if (startX) {
        B1y <- bh_reject(pvals(bk, B0$x), alpha) + dx
        B1x <- bh_reject(pvals(bk, B1y), alpha)
      } else {
        B1x <- bh_reject(pvals(bk, B0$y), alpha)
        B1y <- bh_reject(pvals(bk, B1x), alpha) + dx
      }
      list(x=B1x, y=B1y)
    } else {
      px <- pvals(bk, B0$y) # size dx vector
      py <- pvals(bk, B0$x) # size dy vector

      #Note that the positions returned by bh_reject are already the global
      #numbers.
      return(split(bh_reject(c(px, py), alpha)))
    }
  }

  extract <- function(indx) {
    # Start the extraction from indx
    #
    # First do the initialization at indx, and then repeateadly apply update till a fixed point is found, or one of the sets repeates (forming a non-trivial cycle).
    # In case there is a cycle and if none of the consecutive sets in the cycle are sufficiently disjoint (jaccard >= 0.5), the take the union of all those sets and start the extraction from the conglomerate set. Othewise we say is break is found, and the extraction is terminated (without any result).
    #
    #@return The return value is the extract_res field of the final method results.

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

    B0 <- initialize(indx, pvals_init)

    # Check for dud
    if (length(B0$x)*length(B0$y) <= 1) {
      Dud_count <<- Dud_count + 1
      cat("--initialize was a dud\n")
      return(list("indx" = indx, "report" = "dud"))
    }

    cat(paste0("--obtained initial set of size (", length(B0$x), ", ", length(B0$y), ")\n"))

    # Initializing extraction loop

    chain <- list(B0) #The collection of all sets visited during the extraction
    B1 <- split(integer(0)) #Initial value for B1

    did_it_cycle <- FALSE
    break_or_collapsed <- FALSE
    itCount <- 0

    diagnostics("ExtracionLoopBegins")

    # Extraction loop
    repeat {
      itCount <- itCount + 1
      cat(paste0("--Performing update ", itCount, "...\n"))

      B1 <- update(B0, startX = startsAtX)

      if (length(B1$y) * length(B1$x) == 0) {
        break_or_collapsed <- TRUE
        diagnostics("Collapsed")
        cat("----collapse on at least one side\n")
        break
      }

      jaccards <- rlist::list.mapv(chain, jac(., B1))

      diagnostics("AfterUpdate")

      #End loop if fixed point found.
      if (jac(B1, B0) == 0) break

      # Checking for cycles (4.4.1 in CCME paper)
      if(sum(jaccards == 0) > 0) {
          diagnostics("FoundCycle")

          did_it_cycle <- TRUE
          cat("----cycle found")

          Start <- max(which(jaccards == 0))
          cycle_chain <- chain[Start:length(chain)]

          # Checking for cycle break (4.4.1a)
          seq_pair_jaccards <- rep(0, length(cycle_chain) - 1)
          for (j in seq_along(seq_pair_jaccards)) {
            seq_pair_jaccards[j] <- jac(chain[[Start + j - 1]], chain[[Start + j]])
          }
          if (sum(seq_pair_jaccards > 1 - break_thresh) > 0) {# then break needed
            cat("------break found\n")
            diagnostics("FoundBreak")
            break_or_collapsed <- TRUE
            break
          }

          # Create conglomerate set (and check, 4.4.1b)
          B_J <- conglomerate(cycle_chain)

          B_J_check <- rlist::list.mapv(chain, jaccard(B_J, .))
          if (sum(B_J_check == 0) > 0) {
            cat("------old cycle\n")
            break_or_collapsed <- TRUE
            break
          } else {
            cat("------new cycle\n")
            B1 <- B_J
          }
        }
      cat(paste0("----updated to set of size (", length(B1$x), ", ", length(B1$y), ")\n"))
      cat(paste0("----jaccard to previous = ", round(jaccards[length(jaccards)], 3), "\n"))
      chain <- rlist::list.append(chain, B1)
      B0 <- B1
    }

    diagnostic_info <- diagnostics("EndOfExtract")
    # Storing B_new and collecting update info

    if (break_or_collapsed) {
     Dud_count <<- Dud_count + 1
     return(c(list("indx" = indx,
                 "StableComm" = integer(0),
                 "itCount" = itCount,
                 "did_it_cycle" = did_it_cycle,
                 "report" = "break_or_collapse"), diagnostic_info))
    }

    stableComm <- merge(B1)

    if (mask_extracted) {
      mask(bk, B1$x, B1$y)
    }

    clustered <<- union(clustered, stableComm)
    comms <<- rlist::list.append(comms, B1)



    # Checking overlap with previous sets
    if (length(comms) > 1) {
      OL_check <- rlist::list.mapv(comms, jac(., B1))
      if (sum(OL_check < 1 - OL_thres) > 0) {
        OL_count <<- OL_count + 1
      }
    }

    return(c(list("indx" = indx,
                "StableComm" = stableComm,
                "itCount" = itCount, "did_it_cycle" = did_it_cycle,
                "report" = "complete_extraction"), diagnostic_info))

  }

  #-------------------------------------------------------------------------------
  # Extractions

  indices <- c(Xindx, Yindx)

  # Getting node orders. We know what nodes to start with, 
  # but what is a good order to start in?
  
  if (rank_initial_sets) {
    cat("Ranking initial sets\n\n")
    scores <- rlist::list.mapv(indices, {
      if( is.null(start_nodes) || . %in% start_nodes ) {
        B0 <- initialize(., pval_func=pvals_quick)
        if (length(B0$x)*length(B0$y) > 1) {
          score_quick(bk, B0$x, B0$y)
        } else {
          NA
        }
      } else {
        NA
      }
    })
    #Remove the NAs from the ordering
    extractord <- order(scores, decreasing=TRUE, na.last=NA)
  } else {
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
  }
  

  # Extracting
  
  tryCatch({ extract_res <- lapply(extractord, extract) }, 
           interrupt = function(x) {
             message("Interrup recieved. Terminating extract loop")
             stop_extracting <<- TRUE
           })
  stopped_at <-  rlist::list.findi(extract_res, report == "stop_extracting")
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
  nonNullIndxs <- rlist::list.which(final.sets, length(.) > 0)
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
