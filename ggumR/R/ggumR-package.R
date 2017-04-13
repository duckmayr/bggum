#' ggumR
#'
#' ggumR provides R tools for the generalized graded unfolding model.  
#' @name ggumR
#' @docType package
#' @author  JB Duck-Mayr, Patrick Cunha Silva, Luwei Ying
#' @useDynLib ggumR
#' @importFrom Rcpp sourceCpp
#' @importFrom stats dbeta dnorm pbeta qbeta rbeta
NULL
.onUnload <- function (libpath) {
  library.dynam.unload('ggumR', libpath)
}