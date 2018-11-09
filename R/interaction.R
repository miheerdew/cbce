interaction_gui_non_safe <- function(event, env=parent.frame()) {
  switch(event, 
         "Main:Setup" = {

           env$maxit <- get("max_iterations", envir=env)
           env$fixed_points <- 0
           env$duds <- 0
           env$collapsed <- 0
           env$cycled <- 0
           env$overflowed <- 0
           env$tot <- length(get("extractord", envir=env))
           env$start_time <- Sys.time()
           env$timestamp_plot  <- env$timestamp_update <- Sys.time()
           
           plot_func <- function() {
             #Plot the results from the latest non_dud extraction.
             if(is.null(res <- get("res", envir = env)) || ! (res$itCount > 0)) {
               if(is.null(res <- env$old_res)) {
                 return(NULL)
               }
             }
             
             env$old_res <- res
             
             u <- res$update_info
             itCount <- res$itCount
             jacs <- u$consec_jaccards
             
             char_map <- list(default = 1,
                              cycles = 5,
                              fixed_point = 19,
                              collapse = 8,
                              overflow = 13)
             
             #The plotting characters
             chars = rep(char_map$default, itCount)
             chars[u$disjointed] <- char_map$disjoint
             if(res$fixed_point) {
               chars[itCount] <- char_map$fixed_point
             } else if(res$collapsed) {
               chars[itCount] <- char_map$collapse
             } else if(res$overflow) {
               chars[itCount] <- char_map$overflow
             }
             
             chars[rlist::list.mapv(u$cycle_info, itCount)] <- char_map$cycles
             
             #The labels give sizes 
             labs = rlist::list.mapv(u$sizes, sprintf("%d/%d", .$x, .$y))
             
             #opar <- par()
             #par(mfrow=c(1,2))
             plot(0:itCount, jacs,
                  ylim=c(0,1),
                  xlim=c(0, env$maxit + 2), xlab="Iterations", ylab="Consec jaccards",
                  pch=c(char_map$default, chars), 
                  main=sprintf("Extraction time %s. log(p-value)=%.2f", 
                                          format(res$extraction_time, digits=2),
                                          res$log.pvalue))
             text(0:itCount, jacs, labels=labs, pos=3, cex=0.5, bty="l")
             legend("topright", 
                    title="Events",
                    legend=names(char_map), pch=as.numeric(char_map))
             #if(!all(is.na(u$pvals[[itCount]]))) { 
             #  hist(u$pvals[[itCount]], main="P-values from the last iteration", xlab="pvals")
             #}
           }
           
           env$update_status <- function() {
             
             if(is.null(res <- get("res", envir = env))) {
               return(NULL)
             }
             
             env$fixed_points <- env$fixed_points + res$fixed_point
             env$duds <- env$duds + (res$collapsed && (res$itCount == 0))
             env$collapsed <- env$collapsed + (res$collapsed && res$itCount > 0)
             env$cycled <- env$cycled + (res$cycle_count > 0)
             env$overflowed <- env$overflowed + res$overflowed
             i <- get("i", envir=env)
             
             if(Sys.time() - env$timestamp_update > 1 || i == 1) {
               env$timestamp_update <- Sys.time()
               # Update status every second
               tcltk::tclvalue(env$statusLabel) <- sprintf(
                 "%d iterations, %d: fixed points, %d: duds, %d: collapsed, %d: cycled, %d: overflowed",
                 i, 
                 env$fixed_points, 
                 env$duds, 
                 env$collapsed, 
                 env$cycled, 
                 env$overflowed)
               tcltk::tcl("update")
               tcltk::tkconfigure(env$pb, value=i)
               if(Sys.time() - env$timestamp_plot > 5 || i == 1) {
                env$timestamp_plot <- Sys.time()
                # Update plot every 5 seconds
                tkrplot::tkrreplot(env$plot)
               }
               return(TRUE)
             } else {
               return(FALSE)
             }
           }
           
           # ----- Setup the GUI -----------
           env$win <- tcltk::tktoplevel()
           tcltk::tktitle(env$win) <- get0("TITLE", envir=globalenv(), ifnotfound = "CBCE Progress")
           
           env$stop <- tcltk::tclVar(0)
           env$browser <- tcltk::tclVar(0)
           
           env$statusLabel <- tcltk::tclVar("The extraction status will be shown here")
           
           env$pb <- tcltk2::tk2progress(env$win, length=700)
           tcltk::tkconfigure(env$pb, value=0, maximum=env$tot)
           
           Stop.but <- tcltk::tkbutton(env$win, text = "Stop",
                              command = function() tcltk::tclvalue(env$stop) <- 1)
           Browser.but <- tcltk::tkbutton(env$win, text="Browser",
                              command = function() tcltk::tclvalue(env$browser) <- 1)
           
           Status.lab  <- tcltk2::tk2label(env$win, textvariable = env$statusLabel) 
           
           
           env$plot <- tkrplot::tkrplot(env$win, fun=plot_func, hscale=2, vscale = 1.5)
           
           tcltk::tkgrid(env$plot, columnspan=2)
           tcltk::tkgrid(Status.lab, columnspan=2)
           tcltk::tkgrid(Stop.but, Browser.but)
           tcltk::tkgrid(env$pb, columnspan=2)
           
           
           tcltk::tkfocus(env$win)
         },
         "Main:NextExtraction"={
           
           if(env$update_status()) {
             if(as.numeric(tcltk::tclvalue(env$browser))) {
               tcltk::tclvalue(env$browser) <- 0
               return("browse")
             } 
             
             if (as.numeric(tcltk::tclvalue(env$stop))) {
               return("stop")
             } 
           }
         },
         "Main:Filtering"={
           tcltk::tclvalue(env$statusLabel) <- "Filtering bimodules. Please wait."
           tcltk::tcl("update")
         },
         "Main:End" = {
           tcltk::tkdestroy(env$win)
         }, 
         warning("Unknown event in interaction_gui ", event))
}

#' Provides a GUI for interaction with cbce. Pass 
#' \code{interaction=interaction_gui} to cbce to use the GUI interface.
#' 
#' The global variable \code{TITLE} constrols the title of the GUI.
#' Within the GUI, the \code{stop} button will stop the procedure at any 
#' point and return partial resutls. 
#' 
#' \emph{Advance:} If \code{cbce} is run from RStudio, then the \code{browser} button
#' will run the \code{browser()} command and will allow the user to 
#' step into the method.
#' 
#' @param event The name of the event of the callback
#' @param env The environment to access the method internals and store
#' persistent state
#' @export
interaction_gui <- function(event, env=parent.frame()) {
  tryCatch({
    interaction_gui_non_safe(event, env)
  }, error=function(e) {
    warning("Error detected by GUI:", e)
  }, warning=function(w) {
    warning("Warning detected by GUI:", w)
  })
}

interaction_none <- function(...) TRUE

#' Provides a command line progress indicator for \code{\link{cbce}}
#' Just pass \code{interaction=interaction_cli} as an argument to cbce.
#' 
#' @param event The name of the event of the callback
#' @param env The environment to access the method internals and store
#' persistent state.
#' @export
interaction_cli <- function(event, env=parent.frame()) {
  switch (event,
    "Main:Setup" = {
      env$pb <- utils::txtProgressBar(max=length(get("extractord", envir=env)),
                               style=3)
    },
    "Main:NextExtraction" = {
      utils::setTxtProgressBar(env$pb, get("i", envir=env))
    },
    "Main:End" = {
      close(env$pb)
    }, 
    "Main:Filtering" = {
      NULL
    }, 
    warning("Unknown event in interaction_cli ", event)
  )
}