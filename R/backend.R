#source("chisq.R")
#source("normal.R")
#source("helper.R")

# backends should:
#   add more fields to backend.default return objects
#   set the two_sided field to its meaningfull value
#   implement their respective pvals


backend <- function(X, Y, ...) {
  UseMethod("backend")
}
  
backend.default <- function(X, Y, calc_full_cor=TRUE, init_method=1){
  p <- list(full_xy_cor = if (calc_full_cor) cor(X,Y) else NULL,
            calc_full_cor = calc_full_cor,
            dx = ncol(X), n = nrow(X),
            X = scale(X), Y = scale(Y),
            two_sided = NA) #This must be set by the backend
  class(p) <- paste0("initmethod", init_method)
  p
}

pvals <- function(p, B, ...) {
  #Given a testing set B, return the p-values for the opposite side
  if (length(B) == 0) return(integer(0))
  UseMethod("pvals", p)
}

cors <- function(p, A){
  #Given a set A (either a subset of X nodes, or Y nodes)
  #return the correlation to the opposite side.
  UseMethod("cors", p) 
}

init <- function(p, indx, alpha, ...) {
  UseMethod("init", p)
}

cors.default <- function(p, A){
  if (p$calc_full_cor) {
    if (min(A) > p$dx) {
      #A is in the Y set
      return(p$full_xy_cor[, A - p$dx, drop = FALSE])
    } else {
      #A is in the X set
      return(t(p$full_xy_cor[A, , drop = FALSE]))
    }    
  } else {
    if (min(A) > p$dx) {
      return(crossprod(p$X, p$Y[, A - p$dx, drop = FALSE])/(p$n - 1))
    } else {
      return(crossprod(p$Y, p$X[, A, drop = FALSE])/(p$n - 1))
    }
  }
}

init.initmethod1 <- function(p, indx, alpha, conserv = TRUE) {
  fischer_tranformed_cor <- atanh(as.vector(cors(p, indx))) * sqrt(p$n - 3)
  if (p$two_sided) {
    pvals <- 2 * pnorm(abs(fischer_tranformed_cor), lower.tail = FALSE)
  } else {
    pvals <- pnorm(fischer_tranformed_cor, lower.tail = FALSE)
  }
  successes <- bh_reject(pvals, alpha)
  return(successes)
}