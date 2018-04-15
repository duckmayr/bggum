#' Item Response Function
#' 
#' Plots response functions given alpha, delta, and tau parameters.
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
#' @param sub An optional subtitle for the resulting plot
#' @param layout_matrix An integer matrix dictating the layout of the plot;
#'   the default is a one-column matrix with one element for each item
#' @param option_names An optional character vector giving names for the items'
#'   options; if NULL, generic names (e.g. "Option 1", "Option 2", etc.)
#'   are used
#' @param color A logical vector of length one; if TRUE, each option's line
#'   is printed in a different color (automatically selected via viridis),
#'   while if FALSE (the default), each line is plotted in black
#' @param color_palette A vector of colors to use if color = TRUE;
#'   the default is an eight color, colorblind friendly palette
#'   from Wong (2011)
#'
#' @references Wong, Bang. 2011. "Points of view: Color blindness." Nature
#'   Methods 8:441.
#' 
#' @export
irf <- function(a, d, t, from = -3, to = 3, by = 0.01, sub = "",
                layout_matrix = matrix(1:length(a), ncol = 1, byrow = TRUE),
                option_names = NULL, color = FALSE,
                color_palette = c("#000000", "#e69f00", "#56b4e9", "#009e73",
                                  "#f0e442", "#0072b2", "#d55e00", "#cc79a7")) {
    if ( length(a) != length(d) | (length(a) > 1 & length(a) != length(t)) ) {
        stop("Please provide a, d, and t of the same length.", call. = FALSE)
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
    if ( color ) {
        colors <- color_palette
    }
    else {
        colors <- rep("black", max(K))
    }
    if ( is.null(option_names) ) {
        option_names <- paste("Option", 1:max(K))
    }
    for ( j in 1:length(a) ) {
        y = sapply(th, function(x) ggumProbability(1, x, a[j], d[j], t[[j]]))
        graphics::plot(x = th, y = y, col = colors[1],
                       type = "l", main = paste("Item Response Function", sub),
                       xlab = expression(theta), ylab = expression(P[ij](k)),
                       xlim = c(from, to), ylim=c(0, 1.2))
        for ( k in 2:K[j] ) {
            graphics::lines(th, sapply(th, function(x) {
                ggumProbability(k, x, a[j], d[j], t[[j]])
            }), lty = k, col = colors[k])
        }
        graphics::legend(x = from, y = 1.2, option_names, lty = 1:max(K),
                         horiz = TRUE, bty = "n", col = colors)
    }
    graphics::layout(1) # resets layout
}
