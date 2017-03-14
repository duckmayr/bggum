#' GGUM Likelihood
#' 
#' \code{likelihood} calculates the likelihood of a response matrix
#' 
#' Given parameters for the individuals and test items, and a response matrix
#' recording each individual's response to each test item, \code{likelihood}
#' the product of the probability corresponding to each entry in the matrix.
#' 
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
#' 
#' @return A matrix probabilities for the elements of the response matrix.
likelihood <- function(thetas, responseMatrix, alphas, deltas, taus){
  return(stableProd(probMatrix(thetas, responseMatrix, alphas, deltas, taus)))
}