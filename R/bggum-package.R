#' bggum
#'
#' bggum provides R tools for the Bayesian estimation of generalized graded
#' unfolding model parameters. Please see the vignette
#' (via \code{vignette("bggum")}) for a full in-depth practical guide to
#' Bayesian estimation of GGUM parameters.
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
