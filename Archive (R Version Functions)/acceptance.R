#' Helper function for \code{GGUM}
#' 
#' \code{acceptanceTheta} computes a acceptance ratio for the paramenters used 
#' in the MCMC algorithm for the GGUM.
#' 
#' @param currentValue The current value of the paramenter of interest
#' @param proposedValue The proposed value for the paramenter of interest
#' @param thetas The vector of current values of the items' delta parameters
#' @param alphas The vector of current values of the items' alpha parameters
#' @param deltas The vector of current values of the items' delta parameters
#' @param taus The list of current values of the items' tau vectors
#' @param responseMatrix A numeric matrix with N (the number of individuals)
#'   rows and n (the number of items) columns; each i,j element of the matrix
#'   gives the option chosen by individual i for item j
#' 
#' @return The acceptance ratio for the paramenters used in the MCMC algorithm for the GGUM.
#' @export
acceptanceTheta <- function(currentValue, proposedValue, responseMatrix, 
                            currentValuesTheta, proposedValuesTheta, alphas, deltas, taus, index = 1){
   # Calculate the log-likelihood for the proposed theta
   logllProposed <- likelihoodRow(thetas = proposedValuesTheta,
                                  responseMatrix, alphas, deltas, taus, index = index)  
   if(is.infinite(logllProposed)) return(0)
   # Calculate the log-likelihood for the current theta
   logllCurrent <- likelihoodRow(thetas = currentValuesTheta,
                                     responseMatrix, alphas, deltas, taus, index = index)
   # Calculate the prior for the current theta
   logpriorCurrent <- log(getPrior(param = 'theta', value = currentValue))
   # Calculate the prior for the proposed theta
   logpriorProposed <- log(getPrior(param = 'theta', value = proposedValue))
   # Calculate the acceptance ratio
   acceptanceR <- exp((logllProposed+logpriorProposed)-(logllCurrent+logpriorCurrent))
   # Return the acceptance ratio for theta
   if(is.infinite(acceptanceR)){ 
      return(1)
   }else if(acceptanceR>1){
      acceptanceR <- 1
   }else{
      return(acceptanceR)
   }
}

#' @export
acceptanceAlpha <- function(currentValue, proposedValue, responseMatrix, 
                            currentValuesAlpha, proposedValuesAlpha, thetas, deltas, taus, index = 1){
   # Calculate the prior for the proposed alpha
   logpriorProposed <- log(getPrior(param = 'alpha', value = proposedValue)+0.00001)
   if(is.nan(logpriorProposed))return(0)
   # Calculate the log-likelihood for the current alpha
   logllCurrent <- likelihoodCol(alphas = currentValuesAlpha,
                                     responseMatrix = responseMatrix, 
                                     thetas = thetas, deltas = deltas, taus = taus, index = index)
   # Calculate the prior for the current alpha
   logpriorCurrent <- log(getPrior(param = 'alpha', value = currentValue)+0.00001)
   # Calculate the log-likelihood for the proposed alpha
   logllProposed <- likelihoodCol(alphas = proposedValuesAlpha,
                                      responseMatrix = responseMatrix, 
                                      thetas = thetas, deltas = deltas, taus = taus, index = index)  

   
   # Calculate the acceptance ratio
   acceptanceR <- exp((logllProposed+logpriorProposed)-(logllCurrent+logpriorCurrent))
   # Return the acceptance ratio for alpha
   if(is.infinite(acceptanceR)){ 
      return(1)
   }else if(acceptanceR>1){
      acceptanceR <- 1
   }else{
      return(acceptanceR)
   }
}

#' @export
acceptanceDelta <- function(currentValue, proposedValue, currentValuesDelta, proposedValuesDelta, 
                            responseMatrix, thetas, alphas, taus, index = 1){
   # Calculate the prior for the proposed delta
   logpriorProposed <- log(getPrior(param = 'delta', value = proposedValue)+0.00001)
   if(is.nan(logpriorProposed)) return(0)
   # Calculate the log-likelihood for the current delta
   logllCurrent <- likelihoodCol(deltas = currentValuesDelta,
                                     responseMatrix = responseMatrix, 
                                     thetas = thetas, alphas = alphas, taus = taus, index = index)
   # Calculate the prior for the current delta
   logpriorCurrent <- log(getPrior(param = 'delta', value = currentValue)+0.00001)
   # Calculate the log-likelihood for the proposed delta
   logllProposed <- likelihoodCol(deltas = proposedValuesDelta,
                                      responseMatrix = responseMatrix, 
                                      thetas = thetas, alphas = alphas, taus = taus, index = index)  
   # Calculate the acceptance ratio
   acceptanceR <- exp((logllProposed+logpriorProposed)-(logllCurrent+logpriorCurrent))
   # Return the acceptance ratio for delta
   if(is.infinite(acceptanceR)){ 
      return(1)
   }else if(acceptanceR>1){
      acceptanceR <- 1
   }else{
      return(acceptanceR)
   }
}

#' @export
acceptanceTau <- function(currentValue, proposedValue, currentValuesTau, 
                          proposedValuesTau, responseMatrix, thetas, alphas, deltas, index = 1){
   # Calculate the prior for the proposed tau
   logpriorProposed <- log(getPrior(param = 'tau', value = proposedValue)+0.00001)
   if(is.nan(logpriorProposed)) return(0)
   # Calculate the log-likelihood for the current tau
   logllCurrent <- likelihoodCol(taus = currentValuesTau,
                                     responseMatrix = responseMatrix, 
                                     thetas = thetas, alphas = alphas, deltas = deltas, index = index)
   # Calculate the prior for the current tau
   logpriorCurrent <- log(getPrior(param = 'tau', value = currentValue)+0.00001)
   # Calculate the log-likelihood for the proposed tau
   logllProposed <- likelihoodCol(taus = proposedValuesTau,
                                      responseMatrix = responseMatrix, 
                                      thetas = thetas, alphas = alphas, deltas = deltas, index = index)  
   # Calculate the acceptance ratio
   acceptanceR <- exp((logllProposed+logpriorProposed)-(logllCurrent+logpriorCurrent))
   # Return the acceptance ratio for tau
   if(is.infinite(acceptanceR)){ 
      return(1)
   }else if(acceptanceR>1){
      acceptanceR <- 1
   }else{
      return(acceptanceR)
   }
}