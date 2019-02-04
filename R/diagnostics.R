#' This is the default diagnositcs method passed to cbce.
#' 
#' It collects data about the internal state of the method 
#' and returns it alonside the result.
#' 
#' Currently, the diagnostics collected are:
#' \enumerate{
#' \item Consecutive_jaccards : The jaccard distance at each stage
#'  \item sizes : The sizes at each stage in the extaction
#'  \item disjointed : Whether disjointed at each point.
#' \item cycle : Whether cycled at each point.
#' \item Extraction time.
#' }
#' 
#' @keywords internal
#' @export
diagnostics <- function(event, e=parent.frame()) {
  address <- strsplit(event, ":")[[1]]
  if(address[1] == "Extract") {
    i <- get("itCount", e)
    j <- i + 1
    
    switch(address[2],
           Setup={
             max.len <- get("max_iterations", e) + 1
             
             e$sizes <- rep(list(NULL), max.len)
             e$consec_jaccards <- rep(NULL, max.len)
             e$disjointed <- NULL
             e$cycle_info <- NULL
             e$start_time <- Sys.time()
           },
           LoopBegins={
             B0 <- get("B0", e)
             e$sizes[[j]] <- list(x=length(B0$x), y=length(B0$y))
             e$consec_jaccards[j] <- 1
           },
           Collapsed={
             B1 <- get("B1", e)
             e$sizes[[j]] <- list(x=length(B1$x), y=length(B1$y))
           },
           SizeExceeded={
             B1 <- get("B1", e)
             e$sizes[[j]] <- list(x=length(B1$x), y=length(B1$y))
           },
           AfterUpdate={
             B1 <- get("B1", e)
             jac <- get("dist_to_prev", e)
             
             e$sizes[[j]] <- list(x=length(B1$x), y=length(B1$y))
             e$consec_jaccards[j] <- jac
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
             update_info <- list("consec_jaccards" = e$consec_jaccards[1:j],
                                 "sizes" = e$sizes[1:j],
                                 "disjointed" = e$disjointed,
                                 "cycle_info" = e$cycle_info)
             
             list("update_info" = update_info, 
                  "extraction_time" = Sys.time() - e$start_time)
           },
           NA)
  }
}

#' Don't collect any diagnostics from the method.
#' 
#' @describeIn diagnositcs
#' 
#' @keywords internal
#' @export
diagnostics_none <- function(event, e=parent.frame()) {
  NULL
}