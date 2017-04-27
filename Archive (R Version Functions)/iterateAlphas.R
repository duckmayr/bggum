#' GGUM iterations for Alphas
#' 
#' \code{iterateAlphas} runs the GGUM MCMC algorithm for alphas
#' 
#' @param thetas A numeric vector of length N; each element of the vector is a
#'   respondent's latent trait parameter
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
#' @return A numeric matrix where each i,j entry is the value for alpha j at
#'   iteration i
#' @export
iterateAlphas <- function(thetas, deltas, taus, times, responseMatrix, K){
  m <- ncol(responseMatrix)
  alphas <- rBeta_ab(m, params=list(shape1=1.5, shape2=1.5, a=.25, b=4))
  resultMatrix <- matrix(0, nrow=times+1, ncol=m)
  resultMatrix[1, ] <- alphas
  for ( t in 1:times ) {
    for ( j in 1:m ) {
      if ( t < 1001 ) {
        proposalSD <- 0.15
      } else {
        proposalSD <- sd(resultMatrix[(t-1000):t, j])
      }
      currentVal <- alphas[j]
      proposal <- proposerAlpha(currentVal, proposalSD)
      llCurrent <- llCol(thetas, responseMatrix[ , j], alphas[j], deltas[j],
                         taus[[j]])
      llProposed <- llCol(thetas, responseMatrix[ , j], proposal, deltas[j],
                          taus[[j]])
      if ( is.infinite(llProposed) ) {
        alphas[j] <- proposal
        next
      }
      logpriorCurrent <- log(getPrior(param = 'alpha', value = currentVal))
      logpriorProposed <- log(getPrior(param = 'alpha', value = proposal))
      ratioAlpha <- exp((llProposed + logpriorProposed)
                        - (llCurrent + logpriorCurrent))
      if ( ratioAlpha > 1 ) {
        ratioAlpha <- 1
      }
      if ( sample(c(TRUE, FALSE), 1, prob=c(ratioAlpha, 1-ratioAlpha)) ) {
        alphas[j] <- proposal
      }
    }
    resultMatrix[t+1, ] <- alphas
  }
  return(resultMatrix)
}
