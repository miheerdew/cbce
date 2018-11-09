# base is a backend which implements pvlas_singleton and cors. The other backends can use that if required.

#' An incomplete backend.
#'
#' Other backends build on this. See and use backend instead.
#'
#' @param X Matrix. The data vector for the X side
#' @param Y Matrix. The data vector for the Y side
#' @param cache.size The cache size storing the correlations. Defaults to 0 (don't store anything).
#' @keywords internal
backend.base <- function(X, Y, cache.size=0) {
  p <- list(dx = ncol(X), dy = ncol(Y), n = nrow(X),
            X = scale(X), Y = scale(Y),
            two_sided = FALSE)
  p <- list2env(p)
  p$bk <- new(CorBackend, X, Y, cache.size%/%2, cache.size%/%2)
  class(p) <- c(class(p), "base")
  p
}


cors.base <- function(p, A){
  p$bk$getCor(A)
}

getTstat.base <- function(p, A) {
  if(p$two_sided) {
    p$bk$getSqTstat(A)
  } else {
    rowSums(p$bk$getCor(A))
  }
}