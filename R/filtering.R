#' Filter and summarize the results
#' 
#' The bimodules are clustered based on their jaccard distances, 
#' and the dendrogram is cut to obtain the same number of bimodules as are obtained 
#' by the effective-count formula. For each of the cut subtree an appropriate bimodule
#' is chosen to represent it.
#' 
#' @param extract_res Pass result$extract_res where result is the returned 
#' by the cbce2 procedure
#' @param show.progress Show a progress bar while computing the Jaccard distances
#' @param plot.dendrogram Plot the dendrogram along with the line it is cut at. 
#' @param count.method The method (1 or 2) to use for the calculating the 
#' effective-number of bimodules
#' @param hclust.method The clustering method to use (passed to hclust)
#' 
#' @return 
#' A data frame, each row of which represents the summary of a filtered bimodule. 
#' The columns of the frame are:
#' \itemize{
#'  \item{index} This is the index of the bimodule in extract_res.
#'  \item{x.size, y.size} These are the sizes of the X, Y set of the bimodules
#'  \item{score} The score assigned to the bimodule based on how strong the
#'  correlation is (vs. expected). 
#'  \item{group.size} This bimodule is a representative for a group of bimodules of this size.
#' }
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
#' 
#' df <- filter_and_summarize(res$extract_res)
#' # or just df <- res$filtered_result.df
#' 
#' # The filtered bimodules:
#' bms <- rlist::list.map(res$extract_res[df$index], bimod) 
#'}
#' 
#'@export
#'@importFrom pipeR "%>>%"
#'@import stats
#'@import graphics
filter_and_summarize <- function(extract_res, 
                                 show.progress=TRUE,
                                 plot.dendrogram=FALSE,
                                 count.method=1,
                                 hclust.method="average") {
  
  rlist::list.takeWhile(extract_res, !is.null(.)) %>>%
    rlist::list.update(index.orig=.i) %>>%
    rlist::list.filter(fixed_point) %>>%
      rlist::list.select(bimod, index.orig, log.pvalue) -> ex.fixed
  
  bms <- rlist::list.map(ex.fixed, bimod)
  
  if(length(bms) <= 1) {
    rlist::list.select(ex.fixed, 
                       x.size=length(bimod$x), 
                       y.size=length(bimod$y),
                       score=-log.pvalue, 
                       group=1,
                       group.size=1,
                       index=index.orig) %>>%
      rlist::list.stack() -> df
    return(df)
  } 
  
  Jac <- jacc_matrix(bms, show.progress=show.progress)
  efnum <- switch(count.method,
                    "1" =  effective.num1(bms),
                    "2" = effective.num2(Jac),
                    stop(sprintf('Unknown count.method: "%s"', count.method))
              )
    
  n.tot <- length(bms)
  n.cut <- ceiling(efnum)
    
  hc <- hclust(as.dist(1-Jac), method=hclust.method) 
  grps <- cutree(hc, k=n.cut)
    
  if(plot.dendrogram) {
      height <- hc$height[n.tot - n.cut]
      plot(hc, labels=FALSE)
      abline(h=height, lty=2)
  }

  #Bms, ex.fixed, df.fixed_all are in same orider.
  rlist::list.select(ex.fixed, 
              x.size=length(bimod$x), 
              y.size=length(bimod$y),
              score=-log.pvalue, 
              group=grps[.i],
              index=.i,
              index.orig) %>>% 
    rlist::list.stack() -> df.fixed_all
  
  #Assign a representative to each group
  df.fixed_all %>>%
    dplyr::group_by(group) %>>% 
      dplyr::summarise(index=index[select_rep(bms[index])], count=n()) -> group_rep

  df <- df.fixed_all[group_rep$index, ] 
  df$index <- df$index.orig
  df$group.size <- group_rep$count
  df$index.orig <- NULL
  df
}


effective.num2 <- function(Jac) {
  sum(1/rowSums(Jac))
}

effective.num1 <- function(bimods, show.progress=FALSE) {
  
  total <- 0
  
  if(show.progress) {
    pb <- utils::txtProgressBar(max=length(bimods))
  }
  
  for(i in seq_along(bimods)) {
    b1 <- bimods[[i]]
    len <- rlist::list.map(b1, length(unique(.)))
    Counts <- matrix(0, nrow=len$x, ncol=len$y)
    
    for(b2 in bimods) {
      match.x <- na.omit(match(b2$x, b1$x))
      match.y <- na.omit(match(b2$y, b1$y))
      
      if(length(match.x)*length(match.y) > 0) {
        Counts[match.x, match.y] <- Counts[match.x, match.y] + 1
      }
    }
    
    expected_overlap <- sum(1/Counts)/(len$x*len$y)
    total <- total + expected_overlap
    
    if(show.progress) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  
  total
  
}

select_rep <- function(bimods) {
  #Returns the index of the chosen representative bimodule
  if(length(bimods) <= 1) {
    return(if(length(bimods) == 1) 1 else NULL)
  }
  
  # Map the genes and SNPs to a uniform range.
  
  all.x <- unique(unlist(rlist::list.map(bimods, .$x)))
  all.y <- unique(unlist(rlist::list.map(bimods, .$y)))
  
  bimods <- rlist::list.map(bimods, 
                     list(x=match(.$x, all.x), y=match(.$y, all.y)))
  
  # Count the number of occurances of each pair
  Count.pairs <- matrix(as.integer(0), nrow=length(all.x), ncol=length(all.y))
  
  for(b in bimods) {
    Count.pairs[b$x, b$y] <- Count.pairs[b$x, b$y] + 1 
  }
  
  counts <- rlist::list.mapv(bimods, {
    sum(Count.pairs[.$x, .$y])
  })
  
  # The bimodule with the largest average count is selected.
  which.max(counts)
}
