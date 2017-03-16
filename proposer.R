#' Helper function for \code{GGUM} 
#' 
#' \code{proposer} generates proposed values for parameters.
#' 
#' @param currentValue The current value of the parameter of interest
#' @param paramType A character vector indicating the parameter type; should be
#'   one of 'alpha', 'delta', 'tau', or 'theta'
#' @param item The item for the parameter of interest
#' @param paramVector The vector or list of values for parameter paramType
#' 
#' @return A ratio used in the MCMC algorithm for the GGUM.
proposer <- function(currentValue, paramType, paramVector){
  if (paramType == 'tau') {
    paramVector <- unlist(paramVector)
  }
  return(rnorm(1, currentValue, sd(paramVector)))
}