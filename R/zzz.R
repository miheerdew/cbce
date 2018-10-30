# Just for exporting RCpp things using roxygen

#' @useDynLib cbce, .registration=TRUE
#' @import methods 
#' @importFrom Rcpp loadModule
#' @import Rcpp
#' 
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("cbce", libpath)
}

#'@export CorBackend
loadModule("cbase", TRUE)
