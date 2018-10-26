# base is a backend which implements pvlas_singleton and cors. The other backends can use that if required.

#' An incomplete backend.
#'
#' Other backends build on this. See and use \link{backend} instead.
#'
#' @param X Matrix. The data vector for the X side
#' @param Y Matrix. The data vector for the Y side
#' @param cache.size The cache size storing the correlations. Defaults to 0 (don't store anything).
backend.base <- function(X, Y, cache.size=0){
  #Thes precomputations are stored to compute pvals_singleton for objects of type base, and store useful settings.
  p <- list(dx = ncol(X), dy = ncol(Y), n = nrow(X),
            X = scale(X), Y = scale(Y),
            two_sided = FALSE,
            normal_vector_pval = function(v) stats::pchisq(sum(v^2), df=length(v), lower.tail = FALSE))
  p <- list2env(p)
  
  
  ## Setup the correlation cache.
  c <- cache.size %/% 2
  
  # Can only cache the following x and y
  p$cache_limit.x <- c %/% (p$dy*NUMERIC_SIZE)
  p$cache_limit.y <- c %/% (p$dx*NUMERIC_SIZE)
  
  # Reserve memory for the cache
  p$cache.x <- matrix(numeric(0), nrow=p$dy, ncol=p$cache_limit.x)
  p$cache.y <- matrix(numeric(0), nrow=p$dx, ncol=p$cache_limit.y)
  
  # Holds the X and Y that are cached.
  # and their time stamps
  # The default value of 0 is fine.
  p$time_stamps.x <- p$cached.x <- integer(p$cache_limit.x) 
  p$time_stamps.y <- p$cached.y <- integer(p$cache_limit.y)
  
  # Just statistics
  p$hit_rate.y <- p$hit_rate.x <- 0
  p$count.y <- p$count.x <- 0
    
  class(p) <- c(class(p), "base")
  p
}


#'@import futile.logger
cors.base <- function(p, A){
  # Correct A to use local numbering.
  if(min(A) > p$dx) {
      #A is from Y
      A <- A - p$dx
      
      if(p$cache_limit.y <= 0) {
        # No cache. Calculate and Return the correlations
        return(crossprod(p$X, p$Y[ ,A])/(p$n - 1))
      }
      
      # See which elelments of A are already cached.
      A.find <- match(A, p$cached.y)
      
      #The indices in A which are new/old
      A.new <- which(is.na(A.find))
      A.found <- which(!is.na(A.find))
      
      found.pos <- A.find[A.found]
      
      p$count.y <- p$count.y + 1
      p$hit_rate.y <- p$hit_rate.y + (length(A.found)/length(A) - p$hit_rate.y)/p$count.y
      flog.debug("Hit Y: %f\n", p$hit_rate.y)
      #Fill the correlation matrix.
      
      R <- matrix(numeric(0), nrow=p$dx, ncol=length(A))
      
      #R[, A.found] <- p$cache.y[, found.pos]
      updateColumnsInPlace(R, A.found, p$cache.y[, found.pos, drop=FALSE])
      
      if(length(A.new) > 0) {
        R[, A.new] <- crossprod(p$X, p$Y[, A[A.new]])/(p$n - 1)
      }
      #Update the time stamps
      p$time_stamps.y <- p$time_stamps.y + 1
      p$time_stamps.y[found.pos] <- 0
      
      # The number of items to replace from the cache
      del.count <- min(p$cache_limit.y - length(A.found), length(A.new))
      
      #The positions for the del.count oldest entries
      del.pos <- which(rank(-p$time_stamps.y, ties.method = "random") <= del.count)
      
      #Save del.count of the new entries 
      save.ind <- A.new[seq_len(del.count)]
      p$cached.y[del.pos] <- A[save.ind]
      p$time_stamps.y[del.pos] <- 0
      
      #Finally update the cache with correlations
      #p$cache.y[, del.pos] <- R[, save.ind]
      updateColumnsInPlace(p$cache.y, del.pos, R[, save.ind, drop=FALSE])
      return(R)
  } else {
      #A is from X
      #No need to correct for local indices
      
      if(p$cache_limit.x <= 0) {
        # No cache. Calculate and Return the correlations
        return(crossprod(p$Y, p$X[, A])/(p$n - 1))
      }
      
      # See which elelments of A are already cached.
      A.find <- match(A, p$cached.x)
      
      #The indices in A which are new/old
      A.new <- which(is.na(A.find))
      A.found <- which(!is.na(A.find))
      
      found.pos <- A.find[A.found]
      
      p$count.x <- p$count.x + 1
      p$hit_rate.x <- p$hit_rate.x + (length(A.found)/length(A) - p$hit_rate.x)/p$count.x
      flog.debug("Hit X: %f\n", p$hit_rate.x)
      #Fill the correlation matrix.
      R <- matrix(numeric(0), nrow=p$dy, ncol=length(A))
      
      #R[, A.found] <- p$cache.x[, found.pos]      
      updateColumnsInPlace(R, A.found, p$cache.x[, found.pos, drop=FALSE])

      if(length(A.new) > 0) {
        R[, A.new] <- crossprod(p$Y, p$X[, A[A.new]])/(p$n - 1)
      }
      #Update the time stamps
      p$time_stamps.x <- p$time_stamps.x + 1
      p$time_stamps.x[found.pos] <- 0
      
      # The number of items to replace from the cache
      del.count <- min(p$cache_limit.x - length(A.found), length(A.new))
      
      #The positions for the del.count oldest entries
      del.pos <- which(rank(-p$time_stamps.x, ties.method = "random") <= del.count)
      
      #Save del.count of the new entries 
      save.ind <- A.new[seq_len(del.count)]
      p$cached.x[del.pos] <- A[save.ind]
      p$time_stamps.x[del.pos] <- 0
      
      #Finally update the cache with correlations
      #p$cache.x[, del.pos] <- R[, save.ind]
      updateColumnsInPlace(p$cache.x, del.pos, R[, save.ind, drop=FALSE])
      return(R)
  }
}

#' @describeIn pvals_singleton implementation for the base class
#' @export
pvals_singleton.base <- function(bk, indx, thresh.alpha) {
  # An easy way to calculate p-values from an indx
  fischer_tranformed_cor <- atanh(as.vector(cors(bk, indx))) * sqrt(bk$n - 3)
  if (bk$two_sided) {
    pvals <- 2 * stats::pnorm(abs(fischer_tranformed_cor), lower.tail = FALSE)
  } else {
    pvals <- stats::pnorm(fischer_tranformed_cor, lower.tail = FALSE)
  }
  return(pvals)
}