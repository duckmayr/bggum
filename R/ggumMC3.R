#' GGUM MC3
#'
#' Metropolis Coupled Markov Chain Monte Carlo Sampling for the GGUM
#'
#' @param data A numeric matrix giving the individuals' responses
#' @param sample_iterations A vector of length one giving the number of
#'   iterations the sampler should complete
#' @param swap_interval The period by which to attempt chain swaps;
#'   e.g. if W = 100, a state swap will be proposed between two adjacent chains
#'   every 100 iterations
#' @param n_temps The number of chains
#' @param temps An optional provision of the temperatures for the chains;
#'   if not provided, each temperature T_t for t > 1 is given by
#'   T_{t-1} * (t + 1), and T_1 = 1.
#'
#' @return A numeric matrix giving the parameter values at each iteration
#'   for the cold chain
#'
#' @seealso \code{\link{ggumProbability}}, \code{\link{ggumMCMC}}
#'
#' @export
ggumMC3 <- function(data, sample_iterations, swap_interval,
                    n_temps = length(temps), temps = NULL) {
    if ( is.null(temps) ) {
        if ( n_temps < 2 ) {
            stop(paste("Please provide a vector of temperatures,",
                       "or set n_temps to a number greater than 1."),
                 call. = FALSE)
        }
        temps <- rep(1.0, n_temps)
        for ( t in 2:n_temps ) {
            temps[t] <- 1.0 / (temps[t-1] * (t + 1))
        }
    }
    return(.ggumMC3(data, sample_iterations, n_temps, swap_interval, temps))
}
