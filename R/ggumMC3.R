#' GGUM MC3
#'
#' Metropolis Coupled Markov Chain Monte Carlo Sampling for the GGUM
#'
#' @param data A numeric matrix giving the individuals' responses
#' @param sample_iterations A vector of length one giving the number of
#'   iterations the sampler should complete (default is 10000)
#' @param burn_iterations A vector of length one giving the number of
#'   iterations to burn in (default is 10000)
#' @param sd_tune_iterations A numeric vector of length one; the number of
#'   iterations to use to tune the proposals before the burn-in period
#'   begins (default is 5000). If 0 is given, the proposals are not tuned.
#' @param temp_tune_iterations A numeric vector of length one; if a temperature
#'   schedule is not provided in the \code{temps} argument and
#'   \code{optimize_temps} = TRUE, \code{temp_tune_iterations} gives the number
#'   of iterations to use to tune each temperature before the burn-in period
#'   begins (default is 5000) -- see \code{\link{tune_temperatures}}
#' @param temp_n_draws A numeric vector of length one; if a temperature
#'   schedule is not provided in the \code{temps} argument and
#'   \code{optimize_temps} = TRUE, \code{temp_n_draws} gives the number
#'   of draws from the temperature finding algorithm to calculate each
#'   temperature (default is 2500) -- see \code{\link{tune_temperatures}}
#' @param swap_interval The period by which to attempt chain swaps;
#'   e.g. if swap_interval = 100, a state swap will be proposed between two
#'   adjacent chains every 100 iterations (default is 1)
#' @param flip_interval (Optional) If given, provides the number of iterations
#'   after which the sign of the thetas and deltas should be changed.
#'   For example, if \code{flip_interval = 1000},
#'   every 1000 iterations the theta and delta parameters will be multiplied
#'   by -1 (a valid parameter value change as discussed in Geyer (1991)).
#' @param n_temps The number of chains; should only be given if \code{temps}
#'   is not specified
#' @param temps (Optional) A numeric vector giving the temperatures;
#'   if not provided and \code{optimize_temps = FALSE}, each temperature T_t
#'   for t > 1 is given by 1 + \code{temp_multiplier} * (t-1), and T_1 = 1,
#'   while if \code{optimize_temps = TRUE}, the temperature schedule is
#'   determined according to an optimal temperature finding algorithm
#'   -- see \code{\link{tune_temperatures}}
#' @param optimize_temps A logical vector of length one; if TRUE and a
#'   temperature schedule is not provided in the \code{temps} argument,
#'   an algorithm is run to determine the optimal temperature schedule
#'   (default is TRUE) -- see \code{\link{tune_temperatures}}
#' @param temp_multiplier A numeric vector of length one; if a temperature
#'   schedule is not provided and \code{optimize_temps = FALSE},
#'   controls the differences between temperatures as described in the
#'   description of the \code{temps} argument (default is 0.1)
#' @param proposal_sds (Optional) A list of length four where is element is a
#'   numeric vector giving standard deviations for the proposals;
#'   the first element should be a numeric vector with a standard deviation
#'   for the proposal for each respondent's theta parameter (the latent trait),
#'   the second a vector with a standard deviation for each item's alpha
#'   (discrimination) parameter, the third a vector with a standard deviation
#'   for each item's delta (location) parameter, and the fourth a vector with
#'   a standard deviation for each item's tau (option threshold) parameters.
#'   If not given, the standard deviations are all set to 1.0 before any
#'   tuning begins.
#' @param theta_init (Optional) Either a numeric vector giving an initial value
#'   for each respondent's theta parameter, or a numeric matrix giving an
#'   initial value for each respondent's theta parameter for each parallel chain;
#'   if not given, the initial values are drawn from the prior distribution
#' @param alpha_init (Optional) Either a numeric vector giving an initial value
#'   for each item's alpha parameter, or a numeric matrix giving an
#'   initial value for each item's alpha parameter for each parallel chain;
#'   if not given, the initial values are drawn from the prior distribution
#' @param delta_init (Optional) Either a numeric vector giving an initial value
#'   for each item's delta parameter, or a numeric matrix giving an
#'   initial value for each item's delta parameter for each parallel chain;
#'   if not given, the initial values are drawn from the prior distribution
#' @param tau_init (Optional) Either a list giving an initial value
#'   for each item's tau vector, or a list of lists giving an
#'   initial value for each item's tau vector for each parallel chain;
#'   if not given, the initial values are drawn from the prior distribution
#' @param theta_prior_params A numeric vector of length two;
#'   the mean and standard deviation of theta parameters' prior distribution
#'   (where the theta parameters have a normal prior; the default is 0 and 1)
#' @param alpha_prior_params A numeric vector of length four;
#'   the two shape parameters and a and b values for alpha parameters' prior
#'   distribution (where the alpha parameters have a four parameter beta prior;
#'   the default is 1.5, 1.5, 0.25, and 4)
#' @param delta_prior_params A numeric vector of length four;
#'   the two shape parameters and a and b values for delta parameters' prior
#'   distribution (where the delta parameters have a four parameter beta prior;
#'   the default is 2, 2, -5, and 5)
#' @param tau_prior_params A numeric vector of length four;
#'   the two shape parameters and a and b values for tau parameters' prior
#'   distribution (where the tau parameters have a four parameter beta prior;
#'   the default is 2, 2, -6, and 6)
#' @param recorded_chain An integer vector of length one giving the chain to
#'   record; the default, 1, is the "cold" chain
#' @param return_sds A logical vector of length one; if TRUE, the proposal
#'   standard deviations are stored in an attribute of the returned object
#'   named "proposal_sds." The default is TRUE.
#' @param return_temps A logical vector of length one; if TRUE, the temperatures
#'   of the parallel chains are stored in an attribute of the returned object
#'   named "proposal_temps." The default is TRUE.
#'
#' @return A numeric matrix giving the parameter values at each iteration
#'   for the cold chain
#'
#' @seealso \code{\link{ggumProbability}}, \code{\link{ggumMCMC}},
#'   \code{\link{tune_temperatures}}
#'
#' @export
ggumMC3 <- function(data, sample_iterations = 10000, burn_iterations = 10000,
                    sd_tune_iterations = 5000, temp_tune_iterations = 5000,
                    temp_n_draws = 2500, swap_interval = 1, flip_interval = NA,
                    n_temps = length(temps), temps = NULL,
                    optimize_temps = TRUE, temp_multiplier = 0.1,
                    proposal_sds = NULL,
                    theta_init = NULL, alpha_init = NULL,
                    delta_init = NULL, tau_init = NULL,
                    theta_prior_params = c(0.0, 1.0),
                    alpha_prior_params = c(1.5, 1.5, 0.25, 4.0),
                    delta_prior_params = c(2.0, 2.0, -5.0, 5.0),
                    tau_prior_params = c(2.0, 2.0, -6.0, 6.0),
                    recorded_chain = 1,
                    return_sds = TRUE, return_temps = TRUE) {
    n <- nrow(data)
    m <- ncol(data)
    K <- integer(m)
    if ( is.na(flip_interval) ) {
        flip_interval <- sample_iterations + 1
    }
    for ( j in 1:m ) {
        K[j] = length(unique(na.omit(data[ , j])))
    }
    if ( is.null(theta_init) ) {
        theta_init <- t(sapply(1:n_temps, function(x) {
            init_thetas(n, theta_prior_params[1], theta_prior_params[2])
        }))
    }
    else if ( is.vector(theta_init) ) {
        theta_init <- matrix(rep(theta_init, n_temps), nrow = n_temps, byrow = TRUE)
    }
    if ( is.null(alpha_init) ) {
        alpha_init <- t(sapply(1:n_temps, function(x) {
            init_alphas(m, alpha_prior_params[1], alpha_prior_params[2],
                        alpha_prior_params[3], alpha_prior_params[4])
        }))
    }
    else if ( is.vector(alpha_init) ) {
        alpha_init <- matrix(rep(alpha_init, n_temps), nrow = n_temps, byrow = TRUE)
    }
    if ( is.null(delta_init) ) {
        delta_init <- t(sapply(1:n_temps, function(x) {
            init_deltas(m, delta_prior_params[1], delta_prior_params[2],
                        delta_prior_params[3], delta_prior_params[4])
        }))
    }
    else if ( is.vector(delta_init) ) {
        delta_init <- matrix(rep(delta_init, n_temps), nrow = n_temps, byrow = TRUE)
    }
    if ( is.null(tau_init) ) {
        tau_init <- lapply(1:n_temps, function(x) {
            init_taus(m, tau_prior_params[1], tau_prior_params[2],
                        tau_prior_params[3], tau_prior_params[4], K)
        })
    }
    else if ( is.atomic(tau_init[[1]]) ) {
        tau_init <- lapply(1:n_temps, function(x) tau_init)
    }
    if ( is.null(proposal_sds) ) {
        if ( sd_tune_iterations > 0 ) {
            proposal_sds <- tune_proposals(data, sd_tune_iterations, K,
                                           theta_init[1,], alpha_init[1,],
                                           delta_init[1,], tau_init[[1]],
                                           theta_prior_params, alpha_prior_params,
                                           delta_prior_params, tau_prior_params)
        }
        else {
            proposal_sds <- list(rep(1.0, n), rep(1.0, m), rep(1.0, m), rep(1.0, m))
        }
    }
    if ( is.null(temps) ) {
        if ( n_temps < 2 ) {
            stop(paste("Please provide a vector of temperatures,",
                       "or set n_temps to a number greater than 1."),
                 call. = FALSE)
        }
        if ( optimize_temps ) {
            temps <- tune_temperatures(data, n_temps, temp_tune_iterations,
                                       temp_n_draws, K, proposal_sds,
                                       theta_prior_params, alpha_prior_params,
                                       delta_prior_params, tau_prior_params)
        }
        else{
            temps <- rep(1.0, n_temps)
            for ( t in 2:n_temps ) {
                temps[t] <- 1.0 / (1 + temp_multiplier*(t-1))
            }
        }
    }
    recorded_chain <- recorded_chain - 1 # C++ is 0-indexed
    result <- .ggumMC3(data, sample_iterations, burn_iterations, n_temps,
                       swap_interval, flip_interval, recorded_chain, temps,
                       theta_init, alpha_init,
                       delta_init, tau_init, n, m, K, proposal_sds,
                       theta_prior_params[1], theta_prior_params[2],
                       alpha_prior_params[1], alpha_prior_params[2],
                       alpha_prior_params[3], alpha_prior_params[4],
                       delta_prior_params[1], delta_prior_params[2],
                       delta_prior_params[3], delta_prior_params[4],
                       tau_prior_params[1], tau_prior_params[2],
                       tau_prior_params[3], tau_prior_params[4])
    colnames(result) <- c(paste0("theta", 1:n),
                          paste0("alpha", 1:m),
                          paste0("delta", 1:m),
                          paste(paste0("tau", rep(1:m, times = K-1)),
                                unlist(c(sapply(K-1, seq_len))),
                                sep = "_"))
    class(result) <- c("ggum", "mcmc", class(result))
    attr(result, "mcpar") <- c(1, sample_iterations, 1)
    attr(result, "proposal_sds") <- "if"(return_sds, proposal_sds, NULL)
    attr(result, "temps") <- "if"(return_temps, temps, NULL)
    return(result)
}
