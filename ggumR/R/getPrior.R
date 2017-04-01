#' Get priors for GGUM parameters
#' 
#' \code{getPrior} gives the probability of drawing value from the prior
#' distribution for param.
#' 
#' @param param A character vector indicating the parameter type; should be one
#'   of 'alpha', 'delta', 'tau', or 'theta'
#' @param value The parameter's value
#' 
#' @return The probability of drawing value from a prior distribution for param.
#' @export
getPrior <- function(param, value){
  alphaDistParams <- list(shape1=1.5, shape2=1.5, a=.25, b=4)
  deltaDistParams <- list(shape1=2, shape2=2, a=-5, b=5)
  tauDistParams <- list(shape1=2, shape2=2, a=-6, b=6)
  return(switch(param, 'alpha'=dBeta_ab(value, params=alphaDistParams),
                'delta'=dBeta_ab(value, params=deltaDistParams),
                'tau'=dBeta_ab(value, params=tauDistParams),
                'theta'=dnorm(value)))
}
