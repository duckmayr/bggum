#' Function for \code{GGUM}
#' 
#' \code{sampler} takes in a matrix of all responses and returns the iteration matrix.
#' 
#' @param thetas The vector of current values of the items' delta parameters
#' @param alphas The vector of current values of the items' alpha parameters
#' @param deltas The vector of current values of the items' delta parameters
#' @param taus The list of current values of the items' tau vectors
#' @param responseMatrix A numeric matrix with N (the number of individuals)
#'   rows and n (the number of items) columns; each i,j element of the matrix
#'   gives the option chosen by individual i for item j
#' 
#' @return the iteration matrix with K (the number of iteration stages) rows and n (the number of items) columns.
#' @export
sampler <- function(responseMatrix, times = 7000, K = 4,
                    startTheta = NULL, startAlpha = NULL, startDelta = NULL, startTaus = NULL){
  if(is.null(startTheta)){
    startTheta <- rnorm(n = nrow(responseMatrix), mean = 0, sd = 1)
  }
  if(is.null(startAlpha)){
    startAlpha <- rep(x = 1, times = ncol(responseMatrix))
  }
  if(is.null(startDelta)){
    startDelta <- seq(from = -2.45, to = 2.45, length.out = ncol(responseMatrix))
  }
  if(is.null(startTaus)){
    tauDistParams <-list(shape1=2, shape2=2, a=-6, b=6)
    startTaus <- as.list(sapply(1:10, function(x){
      c(0, rBeta_ab(3, params=tauDistParams))
    }, simplify=FALSE))
  }
  iterate(responseMatrix = responseMatrix,
          times = times,
          K = K,
          startTheta = startTheta, 
          startAlpha = startAlpha, 
          startDelta = startDelta,
          startTaus = startTaus)
}
