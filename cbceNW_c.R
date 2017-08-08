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
  
  
  source("auxiliary.R", local = TRUE)
  source("pvals.R", local = TRUE)
  source("bh_reject.R", local = TRUE)
  source("initialize.R", local = TRUE)
  source("extractNW_c.R", local = TRUE)
  
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
