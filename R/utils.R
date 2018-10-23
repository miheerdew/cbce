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

filter_overlap <- function(comms, tau, inv.length = FALSE) {
  
  K <- length(comms)
  if (inv.length) {
    scores <- 1 / unlist(lapply(comms, length))
  } else {
    scores <- unlist(lapply(comms, length))
  }
  
  
  jaccard_mat0 <- matrix(0, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      jaccard_mat0[i, j] <- length(intersect(comms[[i]], comms[[j]])) / 
        length(comms[[i]])
    }
  }
  
  jaccard_mat <- jaccard_mat0
  diag(jaccard_mat) <- 0
  max_jacc <- max(jaccard_mat)
  deleted_comms <- integer(0)
  
  while (max_jacc > tau) {
    
    inds <- which(jaccard_mat == max_jacc, arr.ind = TRUE)[1, ]
    
    # keep comm with larger score
    delete_comm <- inds[which.min(c(scores[inds[1]], scores[inds[2]]))]
    jaccard_mat[delete_comm, ] <- 0
    jaccard_mat[, delete_comm] <- 0
    deleted_comms <- c(deleted_comms, delete_comm)
    max_jacc <- max(jaccard_mat)
    
  }
  
  kept_comms <- setdiff(1:K, deleted_comms)
  
  return(list("final_comms" = comms[kept_comms],
              "kept_comms" = kept_comms))
  
}

bh_reject <- function(pvals, alpha, multiple_testing_method = 'BHY') {
  
  #pvals may now have NAs, which should not be selected 
  m <- length(pvals)
  if (m == 0) return(integer(0))
  
  ranks <- rank(pvals, ties.method = "first")
  
  pvals_adj <- switch(multiple_testing_method,
        BHY=sum(1 / c(1:m)) * m * pvals / ranks, 
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
  C1 <- list.map(B1, length(unique(.)))
  C2 <- list.map(B2, length(unique(.)))
  
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
