#' @useDynLib bootUR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @import Rcpp

.onUnload <- function (libpath) {
  library.dynam.unload("bootUR", libpath)
}

