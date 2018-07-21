#' bggum
#'
#' bggum provides R tools for the generalized graded unfolding model.
#' @name bggum
#' @docType package
#' @author  JB Duck-Mayr, Patrick Cunha Silva, Luwei Ying
#' @useDynLib bggum
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit
NULL
.onUnload <- function (libpath) {
  library.dynam.unload('bggum', libpath)
}
