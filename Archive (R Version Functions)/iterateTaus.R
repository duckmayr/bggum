#' GGUM iterations for Taus
#' 
#' \code{iterateTaus} runs the GGUM MCMC algorithm for taus
#' 
#' @param thetas A numeric vector of length N; each element of the vector is a
#'   respondent's latent trait parameter
#' @param alphas A numeric vector of length n; each element of the vector is an
#'   item's discrimination parameter
#' @param deltas A numeric vector of length n; each element of the vector is an
#'   item's location parameter
#' @param times A numeric vector of length 1. 
#' @param responseMatrix A numeric matrix with N rows and n (the number of
#'   items) columns; each i,j element of the matrix gives the option chosen by
#'   individual i for item j
#' @param K A numeric vector of length one giving the number of options for
#'   each item
#' 
#' @return A numeric matrix where each i,j entry is the value for
#'   tau[[((j - 1) %/% K) + 1]][((j - 1) %% K) + 1] at iteration i
#' @export
iterateTaus <- function(thetas, alphas, deltas, times, responseMatrix, K){
  m <- ncol(responseMatrix)
  taus <- as.list(sapply(1:m, function(x){
    c(0, rBeta_ab(K-1, params=list(shape1=2, shape2=2, a=-6, b=6)))
    # c(0, runif(K-1, -2, -1))
  }, simplify=FALSE))
  resultMatrix <- matrix(0, nrow=times+1, ncol=m*K)
  resultMatrix[1, ] <- unlist(taus)
  for ( t in 1:times ) {
    for ( j in 1:m ) {
      for ( k in 2:K ) {
        if ( t < 1001 ) {
          proposalSD <- 0.15
        } else {
          proposalSD <- sd(resultMatrix[(t-1000):t, ((j-1)*K)+k])
        }
        currentVal <- taus[[j]][k]
        proposal <- proposerTau(currentVal, proposalSD)
        proposedTauVec <- taus[[j]]
        proposedTauVec[k] <- proposal
        llCurrent <- llCol(thetas, responseMatrix[ , j], alphas[j], deltas[j],
                           taus[[j]])
        llProposed <- llCol(thetas, responseMatrix[ , j], alphas[j], deltas[j],
                            proposedTauVec)
        if ( is.infinite(llProposed) ) {
          taus[[j]][k] <- proposal
          next
        }
        logpriorCurrent <- log(getPrior(param = 'tau', value = currentVal))
        logpriorProposed <- log(getPrior(param = 'tau', value = proposal))
        ratioTau <- exp((llProposed + logpriorProposed)
                          - (llCurrent + logpriorCurrent))
        if ( ratioTau > 1 ) {
          ratioTau <- 1
        }
        if ( sample(c(TRUE, FALSE), 1, prob=c(ratioTau, 1-ratioTau)) ) {
          taus[[j]][k] <- proposal
        }
      }
    }
    resultMatrix[t+1, ] <- unlist(taus)
  }
  return(resultMatrix)
}
