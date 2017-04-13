#' GGUM iterations for Deltas
#' 
#' \code{iterateDeltas} runs the GGUM MCMC algorithm for deltas
#' 
#' @param thetas A numeric vector of length N; each element of the vector is a
#'   respondent's latent trait parameter
#' @param alphas A numeric vector of length n; each element of the vector is an
#'   item's discrimination parameter
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
#' @return A numeric matrix where each i,j entry is the value for delta j at
#'   iteration i
#' @export
iterateDeltas <- function(thetas, alphas, taus, times, responseMatrix, K){
  m <- ncol(responseMatrix)
  deltas <- rBeta_ab(m, params=list(shape1=2, shape2=2, a=-5, b=5))
  resultMatrix <- matrix(0, nrow=times+1, ncol=m)
  resultMatrix[1, ] <- deltas
  for ( t in 1:times ) {
    for ( j in 1:m ) {
      if ( t < 1001 ) {
        proposalSD <- 0.15
      } else {
        proposalSD <- sd(resultMatrix[(t-1000):t, j])
      }
      currentVal <- deltas[j]
      proposal <- proposerDelta(currentVal, proposalSD)
      llCurrent <- llCol(thetas, responseMatrix[ , j], alphas[j], deltas[j],
                         taus[[j]])
      llProposed <- llCol(thetas, responseMatrix[ , j], alphas[j], proposal,
                          taus[[j]])
      if ( is.infinite(llProposed) ) {
        deltas[j] <- proposal
        next
      }
      logpriorCurrent <- log(getPrior(param = 'delta', value = currentVal))
      logpriorProposed <- log(getPrior(param = 'delta', value = proposal))
      ratioDelta <- exp((llProposed + logpriorProposed)
                        - (llCurrent + logpriorCurrent))
      if ( ratioDelta > 1 ) {
        ratioDelta <- 1
      }
      if ( sample(c(TRUE, FALSE), 1, prob=c(ratioDelta, 1-ratioDelta)) ) {
        deltas[j] <- proposal
      }
    }
    resultMatrix[t+1, ] <- deltas
  }
  return(resultMatrix)
}
