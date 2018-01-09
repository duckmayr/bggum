#' Estimate Modes
#' 
#' Estimates one or more modes of a distribution sample.
#' 
#' @param x A numeric vector, the distribution sample
#' @param n A numeric vector of length one, the number of modes to return
#' 
#' @return A numeric vector of length n giving the approximate modes
#' @export
estModes <- function(x, n=1) {
  # Adapted from https://stackoverflow.com/questions/27418461/
  d <- stats::density(x)
  modeInds <- NULL
  for ( i in 2:(length(d$y)-1) ){
    if ( (d$y[i] > d$y[i-1]) & (d$y[i] > d$y[i+1]) ) {
      modeInds <- c(modeInds, i)
    }
  }
  return(d$x[which(d$y %in% sort(d$y[modeInds], decreasing=TRUE)[1:n])])
}
