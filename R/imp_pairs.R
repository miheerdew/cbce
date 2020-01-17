connectivity_threshold <- function(R, eps=1e-2) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package \"igraph\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  R.abs <- abs(R)
  
  th.ub <- 1
  th.lb <- 0
  
  while(th.ub - th.lb > eps) {
    mid <- (th.ub + th.lb)/2
    
    M <- R.abs > mid
    G <- igraph::graph_from_incidence_matrix(M)
    
    if(igraph::is_connected(G)) {
      th.lb <- mid
    } else {
      th.ub <- mid
    }
  }
  
  th.lb
}


