diagnostics1 <- function(event){
  e <- parent.frame()
  switch(event,
  "ExtractionLoopBegins" = {
    e$consec_jaccards <- NULL
    e$mean_jaccards <- NULL
    e$consec_sizes <- list(c(length(e$B0$x), length(e$B0$y)))
    e$found_cycle <- NULL
    e$found_break <- NULL
    e$pvals <- list()
    e$initial_set <- unlist(e$B0)
    e$start_time <- Sys.time()
  }, 
  "AfterUpdate" = {
    consec_jaccard <- jaccard_pairs(e$B0, e$B1)
    e$consec_jaccards <- c(e$consec_jaccards, consec_jaccard)
    #e$mean_jaccards <- c(e$mean_jaccards, mean(e$jaccards))
    e$consec_sizes <- c(e$consec_sizes, list(c(length(e$B1$x), length(e$B1$y))))
    e$found_cycle <- c(e$found_cycle, FALSE)
    e$found_break <- c(e$found_break, FALSE)
  }, 
  "EndOfExtract" = {
    update_info <- list("mean_jaccards" = e$mean_jaccards, 
                        "consec_jaccards" = e$consec_jaccards,
                        "consec_sizes" = e$consec_sizes,
                        "found_cycle" = e$found_cycle,
                        "found_break" = e$found_break,
                        "pvals" = e$pvals)
    list("update_info" = update_info, 
         "initial_set" = e$initial_set,
         "extraction_time" = Sys.time() - e$start_time
         )
  },
  "FoundCycle" = {
    e$found_cycle[e$itCount] <- TRUE
  },
  "FoundBreak" = {
    e$found_break[e$itCount] <- TRUE 
  }, 
  "UpdatePvalues" = {
    f <- parent.frame(n=2)
    qs <- seq(0, 1, length.out = 1000) 
    f$pvals <- c(f$pvals, list(quantile(c(e$px, e$py), qs)))
  }, NA)
}

#' @export
diagnostics2 <- function(event, e=parent.frame()) {
  address <- strsplit(event, ":")[[1]]
  if(address[1] == "Extract") {
    i <- get("itCount", e)
    
    switch(address[2],
           Setup={
             maxit <- get("max_iterations", e)
             
             e$pvals <- rep(list(NULL), maxit)
             e$sizes <- rep(list(NULL), maxit)
             e$consec_jaccards <- rep(NULL, maxit)
             e$disjointed <- NULL
             e$cycle_info <- NULL
             e$start_time <- Sys.time()
           },
           Collapsed={
             B1 <- get("B1", e)
             e$sizes[[i]] <- list(x=length(B1$x), y=length(B1$y))
           },
           AfterUpdate={
             B1 <- get("B1", e)
             jac <- get("dist_to_prev", e)
             
             e$sizes[[i]] <- list(x=length(B1$x), y=length(B1$y))
             e$consec_jaccards[i] <- jac
           },
           
           Disjoint={
             e$disjointed <- c(e$disjointed, i)
           },
           
           FoundCycle={
             chain <- get("chain", e)
             h <- chain[i]
             
             e$cycle_info <- c(e$cycle_info, list(c(itCount=i, length=which.max(h %in% chain[i-1:1]))))
           },
           End={
             update_info <- list("consec_jaccards" = e$consec_jaccards[1:i],
                                 "sizes" = e$sizes[1:i],
                                 "pvals" = e$pvals[1:i],
                                 "disjointed" = e$disjointed,
                                 "cycle_info" = e$cycle_info)
             
             list("update_info" = update_info, 
                  "extraction_time" = Sys.time() - e$start_time)
           },
           NA)
  }
}

#' @export
diagnostics_none <- function(event, e=parent.frame()) {
  NULL
}