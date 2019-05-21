symdiff <- function(s1, s2) {
  return(union(setdiff(s1, s2), setdiff(s2, s1)))
}

jaccard <- function(s1, s2) {
  #The jaccard distance between vectors s1 and s2
  return(length(symdiff(s1, s2)) / length(union(s1, s2)))
}

disjointness <- function(s1, s2) {
  num <- length(intersect(s1, s2))
  den <- min(length(s1), length(s2))
  1 - num/den
}

disjointness_pairs <- function(B1, B2) {
  (disjointness(B1$x, B2$x) + disjointness(B1$y, B2$y))/2
}

intersect_pairs <- function(B1, B2) {
  list(x=intersect(B1$x, B2$x), y=intersect(B1$y, B2$y))  
}

jaccard_pairs <- function(B1, B2) {
  (jaccard(B1$x, B2$x) + jaccard(B1$y, B2$y))/2
}

bh_reject <- function(pvals, alpha, multiple_testing_method = 'BY') {
  
  #pvals may now have NAs, which should not be selected 
  m <- length(pvals)
  if (m == 0) return(integer(0))
  
  ranks <- rank(pvals, ties.method = "first")
  
  pvals_adj <- switch(multiple_testing_method,
        BY=sum(1 / c(1:m)) * m * pvals / ranks, 
        BH=m * pvals / ranks,
        sqrt_BH= m * pvals / sqrt(ranks),
        sqrt_BHY= sum(1 / c(1:m)) * m * pvals / sqrt(ranks),
        square_BH= pvals / (1/m + (ranks/m)^2),
        Bonferroni=m * pvals, 
      )
  
  if (any(pvals_adj <= alpha, na.rm = TRUE)) {
    candidates <- which(pvals_adj <= alpha)
    thres <- max(pvals[candidates])
    return(which(pvals <= thres))
  } else {
    return(integer(0))
  }
  
}

jaccard_sim <- function(B1, B2) {
  #Extension of Jaccard similarity to bimodules by considering it as a set of pairs.
  C12 <- list(x=length(intersect(B1$x, B2$x)), y=length(intersect(B1$y, B2$y)))
  C1 <- rlist::list.map(B1, length(unique(.)))
  C2 <- rlist::list.map(B2, length(unique(.)))
  
  common_pairs <- C12$x*C12$y
  total_pairs <- C1$x*C1$y + C2$x*C2$y - common_pairs
  
  common_pairs/total_pairs
}

jacc_matrix <- function(bimods, show.progress=FALSE) {
  k <- length(bimods)
  
  if(show.progress) {    
    pb <- utils::txtProgressBar(max=choose(k, 2))
  }
  
  Jac <- matrix(numeric(0), nrow=k, ncol=k)
  
  for(i in seq_along(bimods)) {
    for(j in seq_len(i-1)) {
      Jac[i, j] <- Jac[j, i] <- jaccard_sim(bimods[[i]], bimods[[j]])
      if(show.progress) {
        utils::setTxtProgressBar(pb, choose(i-1, 2) + j)
      }
    }
  }
  diag(Jac) <- 1
  Jac
}

NUMERIC_SIZE <- 8