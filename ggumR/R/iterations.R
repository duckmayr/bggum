#' GGUM interations for Thetas, Alphas and Deltas.
#' 
#' \code{MCMC} runs the MCMC algorithm  for GGUM
#' 
#' @param startvalue A numeric vector of length 1.
#' @param iterations A numeric vector of length 1. 
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
iterations <- function(startvalue, iterations, responseMatrix, thetas, alphas, deltas, taus){
   currentTheta <- storeTheta <- chainTheta <- startvalue
   currentAlpha <- storeAlpha <- chainAlpha <- startvalue
   currentDelta <- storeDelta <- chainDelta <- startvalue
   for (i in 1:iterations){
      # Start the BurnIn period. Proposer has sd = 1
      if(i<1001){
        # Do this for Theta
         proposalTheta <- proposer(chainTheta[i])
         ratioTheta <- acceptanceTheta(currentValue = currentTheta, 
                                  proposedValue = proposalTheta,
                                  responseMatrix = responseMatrix, 
                                  alphas = alphas, 
                                  deltas = deltas, 
                                  taus = taus)
         chainTheta[i+1] <- ifelse(sample(c(TRUE, FALSE), 1, prob = c(ratioTheta, 1-ratioTheta)),
                              proposalTheta, chainTheta[i])
         storeTheta[i] <- currentTheta <- chainTheta[i+1]
        # Do this for Alpha
         proposalAlpha <- proposer(chainAlpha[i])
         ratioAlpha <- acceptanceAlpha(currentValue = currentAlpha, 
                                  proposedValue = proposalAlpha,
                                  responseMatrix = responseMatrix, 
                                  thetas = thetas, 
                                  deltas = deltas, 
                                  taus = taus)
         chainAlpha[i+1] <- ifelse(sample(c(TRUE, FALSE), 1, prob = c(ratioAlpha, 1-ratioAlpha)),
                              proposalAlpha, chainAlpha[i])
         storeAlpha[i] <- currentAlpha <- chainAlpha[i+1]         
      }else if(i<5001){ # Still BurnIn period. Proposer has sd = last 1000 thetas
        # Do this for Theta
         proposalTheta <- proposer(chainTheta[i])
         sdPeriodType <- sd(tail(storeTheta, 1000))
         ratioTheta <- acceptanceTheta(currentValue = currentTheta, 
                                      proposedValue = proposalTheta,
                                      responseMatrix = responseMatrix, 
                                      alphas = alphas, 
                                      deltas = deltas, 
                                      taus = taus)
         chainTheta[i+1] <- ifelse(sample(c(TRUE, FALSE), 1, prob = c(ratioTheta, 1-ratioTheta)),
                                  proposalTheta, chainTheta[i])
         storeTheta[i] <- currentTheta <- chainTheta[i+1]
         # Do this for Alpha
         proposalAlpha <- proposer(chainAlpha[i])
         sdPeriodType <- sd(tail(storeAlpha, 1000))
         ratioAlpha <- acceptanceAlpha(currentValue = currentAlpha, 
                                       proposedValue = proposalAlpha,
                                       responseMatrix = responseMatrix, 
                                       thetas = thetas, 
                                       deltas = deltas, 
                                       taus = taus)
         chainAlpha[i+1] <- ifelse(sample(c(TRUE, FALSE), 1, prob = c(ratioAlpha, 1-ratioAlpha)),
                                   proposalAlpha, chainAlpha[i])
         storeAlpha[i] <- currentAlpha <- chainAlpha[i+1]         
      } else { # End of the BurnIn period.
        # Do this for Theta
         sdPeriodType <- sd(tail(storeTheta, 1000))
         proposalTheta <- proposer(chainTheta[i], sdPeriodType = sdPeriodType)
         ratioTheta <- acceptanceTheta(currentValue = currentTheta, 
                                       proposedValue = proposalTheta,
                                       responseMatrix = responseMatrix, 
                                       alphas = alphas, 
                                       deltas = deltas, 
                                       taus = taus)
         chainTheta[i+1] <- ifelse(sample(c(TRUE, FALSE), 1, prob = c(ratioTheta, 1-ratioTheta)),
                                   proposalTheta, chainTheta[i])
         storeTheta[i] <- currentTheta <- chainTheta[i+1]
         # Do this for Alpha
         sdPeriodType <- sd(tail(storeAlpha, 1000))
         proposalAlpha <- proposer(chainAlpha[i], sdPeriodType = sdPeriodType)
         ratioAlpha <- acceptanceAlpha(currentValue = currentAlpha, 
                                       proposedValue = proposalAlpha,
                                       responseMatrix = responseMatrix, 
                                       thetas = thetas, 
                                       deltas = deltas, 
                                       taus = taus)
         chainAlpha[i+1] <- ifelse(sample(c(TRUE, FALSE), 1, prob = c(ratioAlpha, 1-ratioAlpha)),
                                   proposalAlpha, chainAlpha[i])
         storeAlpha[i] <- currentAlpha <- chainAlpha[i+1] 
      }
   }
   return(cbind(chainTheta, chainAlpha))
}
