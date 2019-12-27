#' Post-process a Posterior Sample
#'
#' Post-process the results of \code{\link{ggumMCMC}} or \code{\link{ggumMC3}}
#' using an artificial identifiability constraint (AIC).
#'
#' Since under the GGUM the probability of a response is the same for any given
#' choice of theta and delta parameters and the negative of that choice; i.e.
#'
#' \deqn{Pr(z | \theta, \alpha, \delta, \tau) = Pr(z | -\theta, \alpha, -\delta, \tau),}
#'
#' if symmetric priors are used, the posterior has a reflective mode.
#' This function transforms a posterior sample by enforcing a constraint
#' that a particular parameter is of a given sign, essentially transforming
#' it into a sample from only one of the reflective modes if a suitable
#' constraint is chosen; using a sufficiently extreme parameter is suggested.
#'
#' Please see the vignette (via \code{vignette("bggum")}) for a full in-depth
#' practical guide to Bayesian estimation of GGUM parameters.
#'
#' @param sample A numeric matrix of posterior draws as returned by
#'   \code{\link{ggumMCMC}} or \code{\link{ggumMC3}}.
#' @param constraint An integer vector of length one giving the column number
#'   of the parameter to constrain, or a character vector of length one giving
#'   the column name for the constraint.
#' @param expected_sign A character vector of length one giving the sign for
#'   the constraint; it should be either "-" if the constrained parameter is
#'   to be negative or "+" if the constrained parameter is to be positive.
#'
#' @return A numeric matrix, the post-processed sample.
#'
#' @seealso \code{\link{ggumMCMC}}, \code{\link{ggumMC3}}
#'
#' @examples
#' ## NOTE: This is a toy example just to demonstrate the function, which uses
#' ## a small dataset and an unreasonably low number of sampling interations.
#' ## For a longer practical guide on Bayesian estimation of GGUM parameters,
#' ## please see the vignette ( via vignette("bggum") ).
#' ## We'll simulate data to use for this example:
#' set.seed(123)
#' sim_data <- ggum_simulation(100, 10, 2)
#' ## Now we can generate posterior draws
#' ## (for the purposes of example, we use 100 iterations,
#' ## though in practice you would use much more)
#' draws <- ggumMC3(data = sim_data$response_matrix, n_temps = 2,
#'                  sd_tune_iterations = 100, temp_tune_iterations = 100,
#'                  temp_n_draws = 50,
#'                  burn_iterations = 100, sample_iterations = 100)
#' ## Then you can post-process the output
#' processed_draws <- post_process(sample = draws,
#'                                 constraint = which.min(sim_data$theta),
#'                                 expected_sign = "-")
#'
#' @export
post_process <- function(sample, constraint, expected_sign) {
    if ( ! "ggum" %in% class(sample) ) {
        stop("Provide output from ggumMCMC() or ggumMC3() as sample.")
    }
    comparison <- switch(expected_sign[1],
                         "-" = get(">"),
                         "+" = get("<"),
                         stop("Provide + or - as expected_sign."))
    n <- sum(grepl("theta", colnames(sample)))
    m <- sum(grepl("alpha", colnames(sample)))
    flipcols <- c(1:n, (n+m+1):(n+2*m))
    fliprows <- which(comparison(sample[ , constraint], 0))
    sample[fliprows, flipcols] <- -sample[fliprows, flipcols]
    return(sample)
}

