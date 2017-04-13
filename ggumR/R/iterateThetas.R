#' GGUM iterations for Thetas
#' 
#' \code{iterateThetas} runs the GGUM MCMC algorithm thetas
#' 
#' @param alphas A numeric vector of length n; each element of the vector is an
#'   item's discrimination parameter
#' @param deltas A numeric vector of length n; each element of the vector is an
#'   item's location parameter
#' @param taus A list of numeric vectors; each list element j is a numeric
#'   vector of threshold parameters for item j's options (where the first
#'   element of the vector should be zero).
#' @param times A numeric vector of length 1. 
#' @param responseMatrix A numeric matrix with N rows and n (the number of
#'   items) columns; each i,j element of the matrix gives the option chosen by
#'   individual i for item j
#' @param K A numeric vector of length one giving the number of options for
#'   each item
#' 
#' @return A numeric matrix where each i,j entry is the value for theta j at
#'   iteration i
#' @export
iterateThetas <- function(alphas, deltas, taus, times, responseMatrix, K){
  n <- nrow(responseMatrix)
  m <- ncol(responseMatrix)
  thetas <- rnorm(n)
  resultMatrix <- matrix(0, nrow=times+1, ncol=n)
  resultMatrix[1, ] <- thetas
  for ( t in 1:times ) {
    for ( i in 1:n ) {
      if ( t < 1001 ) {
        proposalSD <- 0.2
      } else {
        proposalSD <- sd(resultMatrix[(t-1000):t, i])
      }
      currentVal <- thetas[i]
      proposal <- proposerTheta(currentVal, proposalSD)
      llCurrent <- llRow(currentVal, responseMatrix[i, ], alphas, deltas, taus)
      llProposed <- llRow(proposal, responseMatrix[i, ], alphas, deltas, taus)
      if ( is.infinite(llProposed) ) {
        thetas[i] <- proposal
        next
      }
      logpriorCurrent <- log(getPrior(param = 'theta', value = currentVal))
      logpriorProposed <- log(getPrior(param = 'theta', value = proposal))
      ratioTheta <- exp((llProposed + logpriorProposed)
                        - (llCurrent + logpriorCurrent))
      if ( ratioTheta > 1 ) {
        ratioTheta <- 1
      }
      if ( sample(c(TRUE, FALSE), 1, prob=c(ratioTheta, 1-ratioTheta)) ) {
        thetas[i] <- proposal
      }
    }
    resultMatrix[t+1, ] <- thetas
  }
  return(resultMatrix)
}
