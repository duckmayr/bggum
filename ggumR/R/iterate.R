#' GGUM iterations for Thetas, Alphas and Deltas.
#' 
#' \code{MCMC} runs the MCMC algorithm  for GGUM
#' 
#' @param startThetas A numeric vector of length 1.
#' @param startAlpha A numeric vector of length 1.
#' @param startDelta A numeric vector of length 1.
#' @param startTaus A numeric vector of length K.
#' @param times A numeric vector of length 1. 
#' @param thetas A numeric vector of length N (the number of respondents); each
#'   each element of the vector is an individual's latent trait parameter
#' @param responseMatrix A numeric matrix with N rows and n (the number of
#'   items) columns; each i,j element of the matrix gives the option chosen by
#'   individual i for item j
#' @param alphas A numeric vector of length n; each element of the vector is an
#'   item's discrimination parameter
#' @param deltas A numeric vector of length n; each element of the vector is an
#'   item's location parameter
#' @param taus A list of numeric vectors; each list element j is a numeric
#'   vector of threshold parameters for item j's options (where the first
#'   element of the vector should be zero).
#' @export
iterate <- function(startTheta, startAlpha, startDelta, startTaus, times, responseMatrix, K){
   n <- nrow(responseMatrix)
   m <- ncol(responseMatrix)
   resultMatrix <- matrix(0, nrow = times+1, ncol = n+2*m+(K)*m)
   resultMatrix[1,] <- c(startTheta, startAlpha, startDelta, unlist(startTaus))
   for (t in 1:times){
        # Do this for Theta
         for(i in 1:n){
         proposalTheta <- proposerTheta(resultMatrix[t,i])
         proposalThetas <- resultMatrix[t, 1:n]
         proposalThetas[i] <- proposalTheta
         ratioTheta <- acceptanceTheta(currentValue = resultMatrix[t,i], 
                                  proposedValue = proposalThetas[i],
                                  proposedValuesTheta = proposalThetas,
                                  currentValuesTheta = resultMatrix[t, 1:n],
                                  responseMatrix = responseMatrix, 
                                  alphas = resultMatrix[t, (n+1):(m+n)], 
                                  deltas = resultMatrix[t, (n+m+1):(n+2*m)], 
                                  taus =  as.list(sapply(1:m, function(x){
                                     resultMatrix[t, (n+2*m+1+(x-1)*K):(n+2*m+K+(x-1)*K)]
                                  }, simplify=FALSE)), 
                                  index = i)
         resultMatrix[t+1, i] <- ifelse(sample(c(TRUE, FALSE), 1, prob = c(ratioTheta, 1-ratioTheta)),
                              proposalTheta, resultMatrix[t,i])
         }
        # Do this for Alpha
         for(j in (n+1):(n+m)){
         proposalAlpha <- proposerAlpha(resultMatrix[t,j])
         proposalAlphas <- resultMatrix[t, (n+1):(n+m)]
         proposalAlphas[i] <- proposalAlpha
         ratioAlpha <- acceptanceAlpha(currentValue = resultMatrix[t,j], 
                                  proposedValue = proposalAlpha,
                                  responseMatrix = responseMatrix,
                                  proposedValuesAlpha = proposalAlphas,
                                  currentValuesAlpha = resultMatrix[t, (n+1):(n+m)],
                                  thetas = resultMatrix[t, 1:n], 
                                  deltas = resultMatrix[t, (n+m+1):(n+2*m)], 
                                  taus =  as.list(sapply(1:m, function(x){
                                     resultMatrix[t, (n+2*m+1+(x-1)*K):(n+2*m+K+(x-1)*K)]
                                  }, simplify=FALSE)),  index = j-n)
         resultMatrix[t+1,j] <- ifelse(sample(c(TRUE, FALSE), 1, prob = c(ratioAlpha, 1-ratioAlpha)),
                              proposalAlpha, resultMatrix[t,j])
         }
         # Do this for Delta
         for(j in (n+m+1):(n+2*m)){
         proposalDelta <- proposerDelta(resultMatrix[t,j])
         proposalDeltas <- resultMatrix[t, (n+m+1):(2*m+n)]
         proposalDeltas[i] <- proposalDelta
         ratioDelta <- acceptanceDelta(currentValue =  resultMatrix[t,j], 
                                       proposedValue = proposalDelta,
                                       proposedValuesDelta = proposalDeltas,
                                       currentValuesDelta = resultMatrix[t, (n+m+1):(2*m+n)],
                                       responseMatrix = responseMatrix, 
                                       thetas = resultMatrix[t, 1:n], 
                                       alphas = resultMatrix[t, (n+1):(m+n)], 
                                       taus =  as.list(sapply(1:m, function(x){
                                          resultMatrix[t, (n+2*m+1+(x-1)*K):(n+2*m+K+(x-1)*K)]
                                       }, simplify=FALSE)), 
                                       index = j-n-m)
         resultMatrix[t+1,j] <- ifelse(sample(c(TRUE, FALSE), 1, prob = c(ratioDelta, 1-ratioDelta)),
                                   proposalDelta, resultMatrix[t,j])
         }
      # Do this for Tau
      for(j in 1:m){
         taus <- as.list(sapply(1:m, function(x){
            resultMatrix[t, (n+2*m+1+(x-1)*K):(n+2*m+K+(x-1)*K)]
         }, simplify=FALSE))
         for(k in 2:K){
            proposalTau <- proposerTau(resultMatrix[t, j])
            proposalTaus <- taus
            proposalTaus[[j]][k] <- proposalTau
            
            ratioTau <- acceptanceTau(currentValue = taus[[j]][k], 
                                      currentValuesTau = taus,
                                      proposedValue = proposalTau,
                                      proposedValuesTau = proposalTaus,
                                      responseMatrix = responseMatrix, 
                                      alphas = resultMatrix[t, (n+1):(m+n)], 
                                      deltas = resultMatrix[t, (n+m+1):(n+2*m)], 
                                      thetas = resultMatrix[t, 1:n], 
                                      index = j)
            taus[[j]][k] <- ifelse(sample(c(TRUE, FALSE), 1, prob = c(ratioTau, 1-ratioTau)),
                                   proposalTau, taus[[j]][k])
         }

         resultMatrix[t+1, (2*m+n+1+(j-1)*(K)):(2*m+n+(j-1)*(K)+(K))] <- unlist(taus[[j]]) 
      }
   }
   return(resultMatrix)
}

#print(c(length((2*m+n+1+(j-1)*(K)):(2*m+n+(j-1)*(K)+(K))), length(unlist(taus))))
