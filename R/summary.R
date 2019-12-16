#' Summarize Posterior Draws for GGUM Parameters
#'
#' Summarize the results of \code{\link{ggumMCMC}} or \code{\link{ggumMC3}}.
#'
#' This function provides the posterior mean, median, standard deviation,
#' and 0.025 and 0.975 quantiles for GGUM parameters from posterior samples
#' drawn using \code{\link{ggumMCMC}} or \code{\link{ggumMC3}}.
#' Please note that the quantiles are calculated using the type 8 algorithm
#' from Hyndman and Fan (1996), as suggested by Hyndman and Fan (1996), rather
#' than the type 7 algorithm that would be the default from R's
#' \code{quantile()}).
#' Before calling this function, care should be taken to ensure that
#' post-processing has been done if necessary to identify the correct
#' reflective posterior mode, as discussed in the vignette and Duck-Mayr and
#' Montgomery (2019).
#'
#' Please see the vignette (via \code{vignette("bggum")}) for a full in-depth
#' practical guide to Bayesian estimation of GGUM parameters.
#'
#' @param object A numeric matrix of posterior draws as returned by
#'   \code{\link{ggumMCMC}} or \code{\link{ggumMC3}}, or a list of
#'   such matrices.
#' @param ... Arguments to be passed to or from other methods
#' @param combine A logical vector of length one; if \code{TRUE} and
#'   \code{object} is a list of \code{ggum} result objects, the matrices are
#'   combined and a summary of the combined sample is given; if \code{FALSE}
#'   and \code{object} is a list of \code{ggum} result objects, each matrix
#'   will be summarized individually; and if \code{object} is not a list, it
#'   has no effect. The default is \code{TRUE}.
#'
#' @return A list with three elements: estimates (a list of length four;
#'   a numeric vector giving the means of the theta draws, a numeric vector
#'   giving the means of the alpha draws, a numeric vector giving the means
#'   of the delta draws, and a list where the means of the tau draws are
#'   collated into a tau estimate vector for each item), sds (a list of length
#'   four giving the posterior standard deviations for the theta, alpha, delta,
#'   and tau draws), and statistics (a matrix with five columns and one row
#'   for each parameter giving the 0.025 quantile, the 0.5 quantile, the mean,
#'   the 0.975 quantile, and the standard deviation of the posterior draws
#'   for each parameter; please note the quantiles are calculated using the
#'   type 8 algorithm from Hyndman and Fan 1996, as suggested by Hyndman and
#'   Fan 1996, rather than the type 7 algorithm that would be the default
#'   from R's \code{quantile()}).
#'
#'   If \code{object} is a list and \code{combine} is \code{FALSE},
#'   a list of such lists will be returned.
#'
#' @seealso \code{\link{ggumMCMC}}, \code{\link{ggumMC3}}
#'
#' @references Duck-Mayr, JBrandon, and Jacob Montgomery. 2019. “Ends Against
#'   the Middle: Scaling Votes When Ideological Opposites Behave the Same for
#'   Antithetical Reasons.” http://jbduckmayr.com/papers/ggum.pdf.
#' @references Hyndman, R. J. and Fan, Y. 1996. "Sample Quantiles in
#'   Packages." \emph{American Statistician} 50, 361--365.
#'
#' @examples
#' \dontrun{
#' ## We'll simulate data to use for this example:
#' set.seed(123)
#' sim_data <- ggum_simulation(1000, 20, 2)
#' ## Now we can generate posterior draws
#' draws <- ggumMC3(data = sim_data$response_matrix, n_temps = 6)
#' ## Then post-process the output
#' processed_draws <- post_process(sample = draws,
#'                                 constraint = which.min(sim_data$theta),
#'                                 expected_sign = "-")
#' ## And now we can obtain a summary of the posterior
#' posterior_summary <- summary(processed_draws)
#' ## It contains all the parameter estimates
#' str(posterior_summary$estimates)
#' ## As well as the posterior standard deviations
#' str(posterior_summary$sds)
#' ## And a matrix of the mean (estimates), median, standard deviations,
#' ## and 0.025 and 0.975 quantiles
#' head(posterior_summary$statistics)
#' }
#'
#' @name summary.ggum
#' @rdname summary.ggum
#' @export
summary.ggum <- function(object, ...) {
    obj_class <- class(object)
    if ( !("ggum" %in% obj_class & "mcmc" %in% obj_class) ) {
        stop(paste("Please provide an object returned by a",
                   "ggumMCMC() or ggumMC() call."), call. = FALSE)
    }
    statistics <- summarize_matrix(object)
    param_names <- colnames(object)
    rownames(statistics) <- param_names
    colnames(statistics) <- c("Quantile 0.025", "Median", "Mean",
                              "Quantile 0.975", "Posterior SD")
    J <- "Mean"
    theta <- statistics[grepl("theta", param_names), J]
    alpha <- statistics[grepl("alpha", param_names), J]
    delta <- statistics[grepl("delta", param_names), J]
    m <- length(delta)
    tau <- lapply(1:m, function(j) {
        c(0, statistics[grepl(paste0("tau", j, "_"), param_names), J])
    })
    tau_names <- unique(sub("_.+", "", param_names[grepl("tau", param_names)]))
    names(tau) <- tau_names
    estimates <- list(theta = theta, alpha = alpha, delta = delta, tau = tau)
    J <- "Posterior SD"
    theta_sds <- statistics[grepl("theta", param_names), J]
    alpha_sds <- statistics[grepl("alpha", param_names), J]
    delta_sds <- statistics[grepl("delta", param_names), J]
    m <- length(delta)
    tau_sds <- lapply(1:m, function(j) {
        c(0, statistics[grepl(paste0("tau", j, "_"), param_names), J])
    })
    names(tau_sds) <- tau_names
    sds <- list(theta_sds = theta_sds, alpha_sds = alpha_sds,
                delta_sds = delta_sds, tau_sds = tau_sds)
    result <- list(estimates = estimates, sds = sds, statistics = statistics)
    class(result) <- "summary.ggum"
    return(result)
}

#' @name summary.ggum
#' @rdname summary.ggum
#' @export
summary.list <- function(object, ..., combine = TRUE) {
    classes <- lapply(object, class)
    if ( all(grepl("ggum", classes)) ) {
        if ( combine ) {
            chain <- do.call("rbind", object)
            class(chain) <- c("ggum", "mcmc")
            result <- summary.ggum(chain)
        } else {
            result <- lapply(object, summary.ggum)
        }
        return(result)
    } else {
        return(summary.default(object, ...))
    }
}
