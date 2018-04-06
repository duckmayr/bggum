#' ggum
#'
#' ggum provides R tools for the generalized graded unfolding model.  
#' @name ggum
#' @docType package
#' @author  JB Duck-Mayr, Patrick Cunha Silva, Luwei Ying
#' @useDynLib ggum
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit
NULL
.onUnload <- function (libpath) {
  library.dynam.unload('ggum', libpath)
}
