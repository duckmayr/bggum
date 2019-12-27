#' Item Characteristic Curve
#'
#' Plots item characteristic curves given alpha, delta, and tau parameters.
#'
#' @param a A numeric vector of alpha parameters
#' @param d A numeric vector of delta parameters
#' @param t Either a list of numeric vectors for the tau parameters for each
#'   option, or a numeric vector if the IRF for only one item is desired --
#'   note the first element of each vector should be zero
#' @param from A numeric vector of length one,
#'   the lowest theta value to estimate response probabilities for;
#'   default is -3
#' @param to A numeric vector of length one,
#'   the highest theta value to estimate response probabilities for;
#'   default is 3
#' @param by A numeric vector of length one giving the spacing between
#'   theta values; default is 0.01
#' @param layout_matrix An integer matrix dictating the layout of the plot;
#'   the default is a one-column matrix with one element for each item
#' @param main_title A character vector giving the plots' main titles;
#'   default is "Item Characteristic Curve".
#' @param sub An optional character vector of subtitles for the resulting plots,
#'   to be pasted onto the main title (helpful for titling individual plots
#'   when plotting multiple items' ICCs).
#' @param color The color to plot the ICC line in; default is "black"
#' @param plot_responses A logical vector of length one specifying whether to
#'   draw points at the theta estimates of actual responses; default is FALSE
#' @param thetas An optional vector of theta estimates for response drawing;
#'   if \code{plot_responses = TRUE} and \code{thetas} is not provided,
#'   an error will be thrown.
#' @param responses An optional matrix or vector (if the ICC for only one item
#'   is desired) of responses; if \code{plot_responses = TRUE} and
#'   \code{responses} is not provided, an error will be thrown.
#'   NOTE: The lowest response for each item should be 0, not 1.
#' @param response_color The color to plot the response points when
#'   \code{plot_responses = TRUE}; the default is \code{"#0000005f"}.
#'
#' @examples
#' ## We'll simulate data to use for these examples:
#' set.seed(123)
#' sim_data <- ggum_simulation(100, 10, 4)
#' ## You can plot the ICC for one item:
#' icc(sim_data$alpha[1], sim_data$delta[1], sim_data$tau[[1]])
#' ## Or multiple items:
#' icc(sim_data$alpha[1:2], sim_data$delta[1:2], sim_data$tau[1:2], sub = 1:2)
#' ## You can also plot the actual responses over the expected response line:
#' icc(sim_data$alpha[1], sim_data$delta[1], sim_data$tau[[1]],
#'     plot_responses = TRUE, responses = sim_data$response_matrix[ , 1],
#'     thetas = sim_data$theta)
#'
#' @export
icc <- function(a, d, t, from = -3, to = 3, by = 0.01, layout_matrix = 1,
                main_title = "Item Characteristic Curve", sub = "",
                color = "black", plot_responses = FALSE,
                thetas = NULL, responses = NULL, response_color = "#0000005f") {
    if ( length(a) != length(d) | (length(a) > 1 & length(a) != length(t)) ) {
        stop("Please provide a, d, and t of the same length.", call. = FALSE)
    }
    m <- length(a)
    if ( plot_responses & is.null(responses) ) {
        stop("Please provide responses when plot_responses = TRUE.")
    }
    if ( plot_responses & is.null(thetas) ) {
        stop("Please provide theta estimates when plot_responses = TRUE.")
    }
    if ( !is.null(responses) ) {
        if ( is.vector(responses) ) {
            responses <- matrix(responses, ncol = 1)
        }
        if ( length(dim(responses)) != 2 ) {
            # In case some random data structure we can't use is passed
            # (e.g. a list rather than a vector, matrix, or data frame)
            stop(paste("Please provide a vector or rectangular data structure",
                       "for responses."), call. = FALSE)
        }
        if ( ncol(responses) != m ) {
            stop("Please provide responses for all items when ",
                 "plot_responses = TRUE.")
        }
    }
    if ( length(main_title) == 1 ) {
        main_title <- rep(main_title, m)
    }
    margins <- "if"(any(main_title != ""), c(3, 3, 3, 1), c(3, 3, 1, 1)) + 0.1
    opar <- graphics::par(mar = margins)
    on.exit(graphics::par(opar))
    if ( length(sub) == 1 ) {
        sub <- rep(sub, m)
    }
    graphics::layout(layout_matrix)
    th <- seq(from = from, to = to, by = by)
    if ( is.list(t) ) {
        K <- sapply(t, length)
    }
    else {
        K <- length(t)
        t <- list(t)
    }
    option_vector <- 1:max(K)
    for ( j in 1:m ) {
        kk <- K[j] - 1
        y <- sapply(th, function(x) {
            sum(0:kk * ggumProbability(0:kk, x, a[j], d[j], t[[j]]))
        })
        graphics::plot(x = th, y = y, col = color, type = "l",
                       main = paste(main_title[j], sub[j]),
                       xlab = "", ylab = "", yaxt = "n", xaxt = "n",
                       xlim = c(from, to), ylim=c(0, kk))
        graphics::axis(side = 2, tick = FALSE, line = -0.75)
        graphics::axis(side = 1, tick = FALSE, line = -0.75)
        graphics::title(ylab = "Expected Response",
                        xlab = expression(theta),
                        line = 1.5)
        if ( plot_responses ) {
            graphics::points(x = thetas, y = responses[ , j], pch = 19,
                             col = response_color)
        }
    }
    graphics::layout(1) # resets layout
}
