#' Helper function for \code{GGUM} 
#' 
#' \code{proposer} generates proposed values for parameters.
#' 
#' @param currentValue The current value of the parameter of interest
#' @param sdPeriodType The standard-deviation used in the t-scaled distribution
#' to draw the proposer.
#' 
#' @return A proposed value for the paramenters used in the MCMC algorithm for the GGUM.

proposer <- function(currentValue, sdPeriodType = 1){
   return(rScaledT(1, df = 1, mu = currentValue, sigma = sdPeriodType))
}