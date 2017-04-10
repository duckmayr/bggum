#' GGUM interations Only Thetas
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
interationsTheta <- function(startvalue, iterations, responseMatrix, thetas, alphas, deltas, taus){
   currentValue <- storevector <- chain <- startvalue
   for (i in 1:iterations){
      # Start the BurnIn period. Proposer has sd = 1
      if(i<1001){
         proposal <- proposer(chain[i])
         ratio <- acceptanceTheta(currentValue = currentValue, 
                                  proposedValue = proposal,
                                  responseMatrix = responseMatrix, 
                                  alphas = alphas, 
                                  deltas = deltas, 
                                  taus = taus)
         chain[i+1] <- ifelse(sample(c(TRUE, FALSE), 1, prob = c(ratio, 1-ratio)),
                              proposal, chain[i])
         storevector[i] <- currentValue <- chain[i+1]
      }else if(i<5001){ # Still BurnIn period. Proposer has sd = last 1000 thetas
         proposal <- proposer(chain[i])
         sdPeriodType <- sd(tail(storevector, 1000))
         ratio <- acceptanceTheta(currentValue = currentValue, 
                                  proposedValue = proposal,
                                  responseMatrix = responseMatrix, 
                                  alphas = alphas, 
                                  deltas = deltas, 
                                  taus = taus)
         chain[i+1] <- ifelse(sample(c(TRUE, FALSE), 1, prob = c(ratio, 1-ratio)),
                              proposal, chain[i])
         storevector[i] <- currentValue <- chain[i+1] 
      } else { # End of the BurnIn period.
         sdPeriodType <- sd(tail(storevector, 1000))
         proposal <- proposer(chain[i], sdPeriodType = sdPeriodType)
         ratio <- acceptanceTheta(currentValue = currentValue, 
                                  proposedValue = proposal,
                                  responseMatrix = responseMatrix, 
                                  alphas = alphas, 
                                  deltas = deltas, 
                                  taus = taus)
         chain[i+1] <- ifelse(sample(c(TRUE, FALSE), 1, prob = c(ratio, 1-ratio)),
                              proposal, chain[i])
         storevector[i] <- currentValue <- chain[i+1]
      }
   }
   return(chain)
}
