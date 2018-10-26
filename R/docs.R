# Just for exporting things using roxygen

#' @useDynLib cbce, .registration=TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("cbce", libpath)
}