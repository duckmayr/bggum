#' GGUM Probability of Response Matrix
#' 
#' \code{probMatrix} gives the probability of each observed response
#' 
#' Given parameters for the individuals and test items, and a response matrix
#' recording each individual's response to each test item, \code{probMatrix}
#' returns a matrix of probabilities; the i,j entry in the matrix is the
#' probability individual i chose the particular response to item j that we
#' observed for individual i and item j (i.e. the probability corresponding to
#' the i,j entry in the response matrix).
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
probMatrix <- function(thetas, responseMatrix, alphas, deltas, taus){
  return(t(sapply(1:length(thetas), function(x){ # for each person
    sapply(1:length(alphas), function(y){ # for each item
      probability(responseMatrix[x, y], length(taus[[y]]), thetas[x],
                  alphas[y],  deltas[y], taus[[y]])})})))
}