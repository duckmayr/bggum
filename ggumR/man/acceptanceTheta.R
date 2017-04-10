#' Helper function for \code{GGUM}
#' 
#' \code{acceptanceTheta} computes a acceptance ratio for Theta used 
#' in the MCMC algorithm for the GGUM.
#' 
#' @param currentValue The current value of theta
#' @param proposedValue The proposed value for theta
#' @param alphas The vector of current values of the items' alpha parameters
#' @param deltas The vector of current values of the items' delta parameters
#' @param taus The list of current values of the items' tau vectors
#' @param responseMatrix A numeric matrix with N (the number of individuals)
#'   rows and n (the number of items) columns; each i,j element of the matrix
#'   gives the option chosen by individual i for item j
#' 
#' @return The acceptance ratio for Theta used in the MCMC algorithm for the GGUM.

acceptanceTheta <- function(currentValue, proposedValue, responseMatrix, alphas, deltas, taus){
   # Calculate the log-likelihood for the current theta
   logllCurrent <- log(likelihoodRow(thetas = currentValue,
                                     responseMatrix, alphas, deltas, taus))
   # Calculate the prior for the current theta
   logpriorCurrent <- log(getPrior(param = 'theta', value = currentValue))
   # Calculate the log-likelihood for the proposed theta
   logllProposed <- log(likelihoodRow(thetas = proposedValue,
                                      responseMatrix, alphas, deltas, taus))      
   # Calculate the prior for the proposed theta
   logpriorCurrent <- log(getPrior(param = 'theta', value = proposedValue))
   # Return the acceptance ratio for theta
   return(exp(logllProposed+logpriorCurrent)-(logllCurrent+logpriorCurrent))
}