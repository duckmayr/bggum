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
#' @param temps (Optional) A numeric vector giving the temperatures;
#'   if not provided, each temperature T_t for t > 1 is given by
#'   T_{t-1} * (t + 1), and T_1 = 1.
#' @param theta_init (Optional) Either a numeric vector giving an initial value
#'   for each respondent's theta parameter, or a numeric matrix giving an
#'   initial value for each respondent's theta parameter for each parallel chain
#' @param alpha_init (Optional) Either a numeric vector giving an initial value
#'   for each item's alpha parameter, or a numeric matrix giving an
#'   initial value for each item's alpha parameter for each parallel chain
#' @param delta_init (Optional) Either a numeric vector giving an initial value
#'   for each item's delta parameter, or a numeric matrix giving an
#'   initial value for each item's delta parameter for each parallel chain
#' @param tau_init (Optional) Either a list giving an initial value
#'   for each item's tau vector, or a list of lists giving an
#'   initial value for each item's tau vector for each parallel chain
#'
#' @return A numeric matrix giving the parameter values at each iteration
#'   for the cold chain
#'
#' @seealso \code{\link{ggumProbability}}, \code{\link{ggumMCMC}}
#'
#' @export
ggumMC3 <- function(data, sample_iterations, swap_interval,
                    n_temps = length(temps), temps = NULL,
                    theta_init = NULL, alpha_init = NULL, delta_init = NULL,
                    tau_init = NULL) {
    n <- nrow(data)
    m <- ncol(data)
    K <- integer(m)
    for ( j in 1:m ) {
        K[j] = length(unique(na.omit(data[ , j])))
    }
    if ( is.null(theta_init) ) {
        theta_init <- t(sapply(1:n_temps, function(x) init_thetas(n, 0.0, 1.5)))
    }
    else if ( is.vector(theta_init) ) {
        theta_init <- matrix(theta_init, nrow = n_temps, byrow = TRUE)
    }
    if ( is.null(alpha_init) ) {
        alpha_init <- t(sapply(1:n_temps, function(x) init_alphas(m, 1.5, 1.5, 0.25, 4.0)))
    }
    else if ( is.vector(alpha_init) ) {
        alpha_init <- matrix(alpha_init, nrow = n_temps, byrow = TRUE)
    }
    if ( is.null(delta_init) ) {
        delta_init <- t(sapply(1:n_temps, function(x) init_deltas(m, 2.0, 2.0, -5.0, 5.0)))
    }
    else if ( is.vector(delta_init) ) {
        delta_init <- matrix(delta_init, nrow = n_temps, byrow = TRUE)
    }
    if ( is.null(tau_init) ) {
        tau_init <- lapply(1:n_temps, function(x) init_taus(m, 2.0, 2.0, -6.0, 6.0, K))
    }
    else if ( is.atomic(tau_init[[1]]) ) {
        tau_init <- lapply(1:n_temps, function(x) tau_init)
    }
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
    return(.ggumMC3(data, sample_iterations, n_temps, swap_interval, temps,
                    theta_init, alpha_init, delta_init, tau_init, K, n, m))
}
