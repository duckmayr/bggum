#' Helper function for \code{GGUM} 
#' 
#' \code{proposer} generates proposed values for parameters.
#' 
#' @param currentValue The current value of the parameter of interest
#' @param paramType A character vector indicating the parameter type; should be
#'   one of 'alpha', 'delta', 'tau', or 'theta'
#' @param param In an unlisted vector of one item's alpha, delta, and tau
#'   parameters, the index of the parameter of interest in this iteration
#' @param paramVector The vector or list of values for parameter paramType
#' 
#' @return A ratio used in the MCMC algorithm for the GGUM.
proposer <- function(currentValue, paramType, param, paramVector){
  if (paramType == 'tau') {
    paramVector <- sapply(paramVector, '[[', param-2)
  }
  return(rnorm(1, currentValue, sd(paramVector)))
}