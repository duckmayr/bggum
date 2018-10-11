#' Summarize Posterior Draws for GGUM Parameters
#'
#' Summarize the results of \code{\link{ggumMCMC}} or \code{\link{ggumMC3}}.
#'
#' @param object A numeric matrix of posterior draws as returned by
#'   \code{\link{ggumMCMC}} or \code{\link{ggumMC3}}
#' @param ... Arguments to be passed to or from other methods
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
#' @seealso \code{\link{ggumMCMC}}, \code{\link{ggumMC3}}
#'
#' @references Hyndman, R. J. and Fan, Y. 1996. "Sample Quantiles in
#'   Packages." American Statistician 50, 361--365.
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
    names(tau) <- sub("_.+", "", param_names[grepl("tau", param_names)])
    estimates <- list(theta = theta, alpha = alpha, delta = delta, tau = tau)
    J <- "Posterior SD"
    theta_sds <- statistics[grepl("theta", param_names), J]
    alpha_sds <- statistics[grepl("alpha", param_names), J]
    delta_sds <- statistics[grepl("delta", param_names), J]
    m <- length(delta)
    tau_sds <- lapply(1:m, function(j) {
        c(0, statistics[grepl(paste0("tau", j, "_"), param_names), J])
    })
    names(tau_sds) <- sub("_.+", "", param_names[grepl("tau", param_names)])
    sds <- list(theta_sds = theta_sds, alpha_sds = alpha_sds,
                delta_sds = delta_sds, tau_sds = tau_sds)
    result <- list(estimates = estimates, sds = sds, statistics = statistics)
    class(result) <- "summary.ggum"
    return(result)
}
