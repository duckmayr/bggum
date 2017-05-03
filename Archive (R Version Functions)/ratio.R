#' Helper function for \code{GGUM}
#' 
#' \code{ratio} computes a ratio used in the MCMC algorithm for the GGUM.
#' 
#' @param param In an unlisted vector of one item's alpha, delta, and tau
#'   parameters, the index of the parameter of interest in this iteration
#' @param person The index or indices of the individuals whose response
#'   probabilities would be affected by accepting the proposed change
#' @param item The index or indices of the items for which the probability of
#'   individuals' responses would be affected by accepting the proposed change
#' @param currentValue The current value of the parameter of interest
#' @param proposedValue The proposed value for the parameter of interest
#' @param alphas The vector of current values of the items' alpha parameters
#' @param deltas The vector of current values of the items' delta parameters
#' @param taus The list of current values of the items' tau vectors
#' @param thetas The vector of current values of the theta parameters
#' @param responseMatrix A numeric matrix with N (the number of individuals)
#'   rows and n (the number of items) columns; each i,j element of the matrix
#'   gives the option chosen by individual i for item j
#' @param paramType A character vector indicating the parameter type; should be
#'   one of 'alpha', 'delta', 'tau', or 'theta'
#' 
#' @return A ratio used in the MCMC algorithm for the GGUM.
ratio <- function(param, person, item, currentValue, proposedValue, alphas,
                  deltas, taus, thetas, responseMatrix, paramType){
  # First we calculate probabilities for the response matrix given the current
  # values of the parameters:
  probMat <- probMatrix(thetas, responseMatrix, alphas, deltas, taus)
  # And the likelihood of the response matrix:
  currentLikelihood <- stableProd(probMat)
  # Then we will replace the entries of probMat at issue given the parameter of
  # interest in this iteration:
  if (paramType == 'theta') {
    probMat[person, ] <- sapply(1:ncol(probMat), function(x){
      probability(responseMatrix[person, x], length(taus[[x]]), proposedValue,
                  alphas[x], deltas[x], taus[[x]])
    })
  } else {
    params <- unlist(c(alphas[item], deltas[item], taus[[item]]))
    params[param] <- proposedValue
    probMat[ , item] <- sapply(1:nrow(probMat), function(x){
      probability(responseMatrix[x, item], length(taus[[item]]), thetas[x],
                  params[1], params[2], params[3:length(params)])
    })
  }
  # And get the likelihood of the response matrix given the proposed change:
  proposedLikelihood <- stableProd(probMat)
  # Then we can calculate the ratio used in the MCMC algorithm:
  return((proposedLikelihood*getPrior(paramType, proposedValue))
         /(currentLikelihood*getPrior(paramType, currentValue)))
}