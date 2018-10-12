plot.progress <- function(...)	{
  vectOfBar <- c(...)*100
  labels <- names(vectOfBar)
  numOfBar <- length(vectOfBar)
  plot(c(0,100), c(0,numOfBar), type='n', xlab='', ylab='', yaxt='n', mar=c(3,3,3,3))
  for(i in 1:numOfBar) {
    rect(0, 0.1+i-1, vectOfBar[i], 0.9+i-1, col=rainbow(numOfBar)[i])
    text(0.5, 0.5+i-1, paste0(labels[i], ': ', round(vectOfBar[i],2), '%'), adj=0)
  }
  title('Progress...')
}

#' @import tcltk2
#' @import tcltk
#' @export
interaction_gui <- function(event, env=parent.frame()) {
  switch(event, 
         "Main:Setup" = {
           
           env$maxit <- get("max_iterations", envir=env)
           env$fixed_points <- 0
           env$duds <- 0
           env$collapsed <- 0
           env$cycled <- 0
           env$overflowed <- 0
           env$tot <- length(get("extractord", envir=env))
           env$timestamp <- env$start_time <- Sys.time()
           
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
             if(length(labs) == 0) {
               msg <- sprintf("Size list: 0, but itCount: %d. Not plotting\n", itCount)
               return(plot(0,0, main = msg))
             } 
             
             #opar <- par()
             #par(mfrow=c(1,2))
             plot(1:itCount, jacs,
                  ylim=c(0,1),
                  xlim=c(0, env$maxit + 2), xlab="Iterations", ylab="Consec jaccards",
                  pch=chars, main=sprintf("Extraction time %s. log(p-value)=%.2f", 
                                          format(res$extraction_time, digits=2),
                                          res$log.pvalue))
             text(1:itCount, jacs, labels=labs, pos=3, cex=0.5, bty="l")
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
             
             if(Sys.time() - env$timestamp > 1) {
               env$timestamp <- Sys.time()
               # Update status every second
               tclvalue(env$statusLabel) <- sprintf(
                 "%d iterations, %d: fixed points, %d: duds, %d: collapsed, %d: cycled, %d: overflowed",
                 i, 
                 env$fixed_points, 
                 env$duds, 
                 env$collapsed, 
                 env$cycled, 
                 env$overflowed)
               tcl("update")
               
               tkconfigure(env$pb, value=i)
               tkrplot::tkrreplot(env$plot)
               return(TRUE)
             } else {
               return(FALSE)
             }
           }
           
           # ----- Setup the GUI -----------
           env$win <- tktoplevel()
           tktitle(env$win) <- TITLE
           
           env$stop <- tclVar(0)
           env$browser <- tclVar(0)
           
           env$statusLabel <- tclVar("The extraction status will be shown here")
           
           env$pb <- tk2progress(env$win, length=700)
           tkconfigure(env$pb, value=0, maximum=env$tot)
           
           Stop.but <- tkbutton(env$win, text = "Stop",
                                command = function() tclvalue(env$stop) <- 1)
           Browser.but <- tkbutton(env$win, text="Browser",
                                   command = function() tclvalue(env$browser) <- 1)
           
           Status.lab  <- tk2label(env$win, textvariable = env$statusLabel) 
           
           
           env$plot <- tkrplot::tkrplot(env$win, fun=plot_func, hscale=2, vscale = 1.5)
           
           tkgrid(env$plot, columnspan=2)
           tkgrid(Status.lab, columnspan=2)
           tkgrid(Stop.but, Browser.but)
           tkgrid(env$pb, columnspan=2)
           
           
           tkfocus(env$win)
         },
         "Main:NextExtraction"={
           
           if(env$update_status()) {
             if(as.numeric(tclvalue(env$browser))) {
               tclvalue(env$browser) <- 0
               return("browse")
             } 
             
             if (as.numeric(tclvalue(env$stop))) {
               return("stop")
             } 
           }
         },
         
         "Main:End" = {
           tkdestroy(env$win)
         }, NA)
}

#' @export
interaction_gui_safe <- function(...) {
  tryCatch({
    interaction_gui(...)
  }, error=function(e) {
    cat(paste("Error:", e))
  }, warning=function(w) {
    cat(paste("Error", e))
  })
}

#' @export
interaction_none <- function(...) TRUE