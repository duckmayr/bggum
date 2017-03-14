#' Helper function for \code{GGUM}
#' 
#' \code{decider} picks one parameter value in a GGUM MCMC algorithm iteration.
#' 
#' @param param In an unlisted vector of one item's alpha, delta, and tau
#'   parameters, the index of the parameter of interest in this iteration
#' @param person The index or indices of the individuals whose response
#'   probabilities would be affected by accepting the proposed change
#' @param item The index or indices of the items for which the probability of
#'   individuals' responses would be affected by accepting the proposed change
#' @param currentValue The current value of the parameter of interest
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
decider <- function(param, person, item, currentValue, alphas, deltas, taus,
                      thetas, responseMatrix, paramType){
  proposedValue <- proposer(currentValue, paramType, param,
                            switch(paramType, 'theta'=thetas, 'alpha'=alphas,
                                   'delta'=deltas, 'tau'=taus))
  reSetProb <- ratio(param, person, item, currentValue, proposedValue, alphas,
                     deltas, taus, thetas, responseMatrix, paramType)
  if (reSetProb > 1) {reSetProb <- 1}
  reSet <- sample(c(TRUE, FALSE), 1, prob=c(reSetProb, 1-reSetProb))
  return(ifelse(reSet, proposedValue, currentValue))
}