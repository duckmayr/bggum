#' GGUM Simulation
#'
#' Generates randomly drawn item and person parameters, and simluated responses.
#'
#' @param n An integer vector of length one giving the number of respondents
#' @param m An integer vector of length one giving the number of items
#' @param K An integer vector giving the number of options for each item;
#'   if the vector is of length one, all m items will have the same number of
#'   options.
#' @param theta (Optional) A numeric vector of respondents' latent traits;
#'   if not given, the values are drawn from a normal distribution whose
#'   mean and standard deviation are given by the \code{theta_params} parameter
#' @param alpha (Optional) A numeric vector of items' discrimination parameters;
#'   if not given, the values are drawn from a four parameter beta distribution
#'   whose parameters are given by \code{alpha_params}
#' @param delta (Optional) A numeric vector of items' location parameters;
#'   if not given, the values are drawn from a four parameter beta distribution
#'   whose parameters are given by \code{delta_params}
#' @param tau (Optional) A list of numeric vectors giving each item's option
#'   thresholds; if not given, the values are drawn from a four parameter beta
#'   distribution whose parameters are given by \code{tau_params}
#' @param theta_params A numeric vector of length two; the mean and standard
#'   deviation of the normal distribution theta is drawn from
#' @param alpha_params A numeric vector of length four; the two shape
#'   parameters and a and b values for the four parameter beta distribution
#'   alpha is drawn from; the default is 1.5, 1.5, 0.25, and 4
#' @param delta_params A numeric vector of length four; the two shape
#'   parameters and a and b values for the four parameter beta distribution
#'   delta is drawn from; the default is 2, 2, -5, and 5
#' @param tau_params A numeric vector of length four;the two shape
#'   parameters and a and b values for the four parameter beta distribution
#'   each tau vector is drawn from; the default is 1.5, 1.5, -2, and 0
#'
#' @return A list with five elements; "theta" containing the theta draws,
#'   "alpha" containing the alpha draws, "delta" containing the delta draws,
#'   "tau" containing the tau draws, and "resp_mat" containing the simulated
#'   response matrix.
#'
#' @seealso \code{\link{ggumProbability}}
#' @export
ggum_simulation <- function(n, m, K, theta = NULL, alpha = NULL, delta = NULL,
                     tau = NULL, theta_params = c(0.0, 1.0),
                     alpha_params = c(1.5, 1.5, 0.25, 4.0),
                     delta_params = c(2.0, 2.0, -5.0, 5.0),
                     tau_params = c(1.5, 1.5, -2.0, 0.0)) {
    if ( length(K) == 1 ) {
        K <- rep(K, m)
    }
    if ( is.null(theta) ) {
        theta <- init_thetas(n, theta_params[1], theta_params[2])
    }
    if ( is.null(alpha) ) {
        alpha <- init_alphas(m, alpha_params[1], alpha_params[2],
                             alpha_params[3], alpha_params[4])
    }
    if ( is.null(delta) ) {
        delta <- init_deltas(m, delta_params[1], delta_params[2],
                             delta_params[3], delta_params[4])
    }
    if ( is.null(tau) ) {
        tau <- init_taus(m, tau_params[1], tau_params[2],
                         tau_params[3], tau_params[4], K)
    }
    resp_mat <- .ggum_simulation(n, m, K, theta, alpha, delta, tau)
    return(list(theta = theta, alpha = alpha, delta = delta, tau = tau,
                resp_mat = resp_mat))
}
