symdiff <- function(s1, s2) {
  return(union(setdiff(s1, s2), setdiff(s2, s1)))
}

jaccard <- function(s1, s2) {
  return(length(symdiff(s1, s2)) / length(union(s1, s2)))
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

bh_reject <- function(pvals, alpha, conserv = TRUE) {
  
  m <- length(pvals)
  if (m == 0) return(integer(0))
  
  if (!conserv) {
    pvals_adj <- m * pvals / rank(pvals, ties.method = "first")
  } else {
    mults <- sum(1 / c(1:m))
    pvals_adj <- mults * m * pvals / rank(pvals, ties.method = "first")
  }
  
  if (sum(pvals_adj <= alpha) > 0) {
    thres <- max(pvals[pvals_adj <= alpha])
    return(which(pvals <= thres))
  } else {
    return(integer(0))
  }
  
}

