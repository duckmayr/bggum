#' bggum
#'
#' bggum provides R tools for the Bayesian estimation of generalized graded
#' unfolding model (Roberts, Donoghue, and Laughlin 2000) parameters.
#' Please see the vignette (via \code{vignette("bggum")}) for a practical
#' guide to Bayesian estimation of GGUM parameters.
#' Duck-Mayr and Montgomery (2019) provides a more detailed theoretical
#' discussion of Bayesian estimation of GGUM parameters.
#'
#' @name bggum
#' @docType package
#' @author  JBrandon Duck-Mayr and Jacob Montgomery
#'
#' @references Duck-Mayr, JBrandon, and Jacob Montgomery. 2019.
#'   \dQuote{Ends Against the Middle: Scaling Votes When Ideological Opposites
#'   Behave the Same for Antithetical Reasons.}
#'   \url{http://jbduckmayr.com/papers/ggum.pdf}.
#' @references Roberts, James S., John R. Donoghue, and James E. Laughlin. 2000.
#'   \dQuote{A General Item Response Theory Model for Unfolding Unidimensional
#'   Polytomous Responses.} \emph{Applied Psychological Measurement}
#'   24(1): 3--32.
#' 
#' @useDynLib bggum
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit
NULL
.onUnload <- function (libpath) {
    library.dynam.unload('bggum', libpath)
}
