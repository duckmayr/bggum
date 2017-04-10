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
   # Calculate the acceptance ratio
   acceptanceR <- exp(logllProposed+logpriorCurrent)-(logllCurrent+logpriorCurrent)
   # Return the acceptance ratio for alpha
   if(is.infinite(acceptanceR)){
      return(1)
   }else {
      return(acceptaceR)
   }
}

#' @export
acceptanceAlpha <- function(currentValue, proposedValue, responseMatrix, thetas, deltas, taus){
   # Calculate the log-likelihood for the current alpha
   logllCurrent <- log(likelihoodCol(alphas = currentValue,
                                     responseMatrix = responseMatrix, 
                                     thetas = thetas, deltas = deltas, taus = taus))
   # Calculate the prior for the current alpha
   logpriorCurrent <- log(getPrior(param = 'alpha', value = currentValue))
   # Calculate the log-likelihood for the proposed alpha
   logllProposed <- log(likelihoodCol(alpha = proposedValue,
                                      responseMatrix = responseMatrix, 
                                      thetas = thetas, deltas = deltas, taus = taus))  
   # Calculate the prior for the proposed alpha
   logpriorCurrent <- log(getPrior(param = 'alpha', value = proposedValue))
   
   # Calculate the acceptance ratio
   acceptanceR <- exp(logllProposed+logpriorCurrent)-(logllCurrent+logpriorCurrent)
   # Return the acceptance ratio for alpha
   if(is.infinite(acceptanceR)){ return(1)
   }else{
      return(acceptanceR)
   }
}

#' @export
acceptanceDelta <- function(currentValue, proposedValue, responseMatrix, thetas, alphas, taus){
   # Calculate the log-likelihood for the current delta
   logllCurrent <- log(likelihoodCol(deltas = currentValue,
                                     responseMatrix = responseMatrix, 
                                     thetas = thetas, alphas = alphas, taus = taus))
   # Calculate the prior for the current delta
   logpriorCurrent <- log(getPrior(param = 'delta', value = currentValue))
   # Calculate the log-likelihood for the proposed delta
   logllProposed <- log(likelihoodCol(deltas = proposedValue,
                                      responseMatrix = responseMatrix, 
                                      thetas = thetas, alphas = alphas, taus = taus))  
   # Calculate the prior for the proposed delta
   logpriorCurrent <- log(getPrior(param = 'delta', value = proposedValue))
   
   # Calculate the acceptance ratio
   acceptanceR <- exp(logllProposed+logpriorCurrent)-(logllCurrent+logpriorCurrent)
   # Return the acceptance ratio for delta
   if(is.infinite(acceptanceR)){ 
      return(1)
   }else{
      return(acceptanceR)
   }
}

#' @export
acceptanceTau <- function(currentValue, proposedValue, responseMatrix, thetas, alphas, deltas){
   # Calculate the log-likelihood for the current tau
   logllCurrent <- log(likelihoodCol(taus = currentValue,
                                     responseMatrix = responseMatrix, 
                                     thetas = thetas, alphas = alphas, deltas = deltas))
   # Calculate the prior for the current tau
   logpriorCurrent <- log(getPrior(param = 'tau', value = currentValue))
   # Calculate the log-likelihood for the proposed tau
   logllProposed <- log(likelihoodCol(taus = proposedValue,
                                      responseMatrix = responseMatrix, 
                                      thetas = thetas, alphas = alphas, deltas = deltas))  
   # Calculate the prior for the proposed tau
   logpriorCurrent <- log(getPrior(param = 'tau', value = proposedValue))
   
   # Calculate the acceptance ratio
   acceptanceR <- exp(logllProposed+logpriorCurrent)-(logllCurrent+logpriorCurrent)
   # Return the acceptance ratio for tau
   if(is.infinite(acceptanceR)){ 
      return(1)
   }else{
      return(acceptanceR)
   }
}
