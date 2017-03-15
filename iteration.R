#' Iteration of the GGUM MCMC algorithm
#' 
#' \code{iteration} performs one iteration of the MCMC algorithm for GGUM
#' 
#' @param thetas The vector of current values of the theta parameters
#' @param alphas The vector of current values of the items' alpha parameters
#' @param deltas The vector of current values of the items' delta parameters
#' @param taus The list of current values of the items' tau vectors
#' @param responseMatrix A numeric matrix with N (the number of individuals)
#'   rows and n (the number of items) columns; each i,j element of the matrix
#'   gives the option chosen by individual i for item j
#' 
#' @return A list of the resulting theta, alpha, delta, and tau values.
iteration <- function(thetas, alphas, deltas, taus, responseMatrix){
  thetas <- sapply(1:length(thetas), function(x){
    decider(0, x, 1:length(alphas), thetas[x], alphas, deltas, taus, thetas,
            responseMatrix, 'theta')
  })
  alphas <- sapply(1:length(alphas), function(x){
    decider(1, 1:length(thetas), x, alphas[x], alphas, deltas, taus, thetas,
            responseMatrix, 'alpha')
  })
  deltas <- sapply(1:length(alphas), function(x){
    decider(2, 1:length(thetas), x, deltas[x], alphas, deltas, taus, thetas,
            responseMatrix, 'delta')
  })
  taus <- sapply(1:length(taus), function(x){
    sapply(2:length(taus[[x]]), function(y){
      decider(y+2, 1:length(thetas), x, taus[[x]][y], alphas, deltas, taus,
              thetas, responseMatrix, 'tau')
    })
  })
  return(list(thetas=thetas, alphas=alphas, deltas=deltas, taus=taus))
}
