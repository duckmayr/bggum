#' bggum
#'
#' bggum provides R tools for the Bayesian estimation of generalized graded
#' unfolding model parameters. For more information, see vignette("bggum").
#'
#' @name bggum
#' @docType package
#' @author  JBrandon Duck-Mayr and Jacob Montgomery
#' @useDynLib bggum
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit
NULL
.onUnload <- function (libpath) {
    library.dynam.unload('bggum', libpath)
}
