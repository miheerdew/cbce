symdiff <- function (s1, s2) {
  return(union(setdiff(s1, s2), setdiff(s2, s1)))
}

jaccard <- function (s1, s2) {
  return(length(symdiff(s1, s2)) / length(union(s1, s2)))
}

filter_overlap <- function (comms, tau, inv.length = FALSE) {
  
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

list2df <- function(listoflists) {
  #Usage: list2df(results$extraction_res)
  field_names <- unique(unlist(lapply(listoflists, names)))
  N <- length(listoflists)
  
  transpose_list <- rep(list(rep(list(), N)), length(field_names))
  names(transpose_list) <- field_names
  
  for(i in 1:N){
    for(n in field_names){
      transpose_list[[n]][i] <- listoflists[[i]][[n]] 
    }
  }
  
  do.call(data.frame, c(lapply(transpose_list, I), list(stringsAsFactors=FALSE)))
}