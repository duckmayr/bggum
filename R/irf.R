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
#' @param layout_matrix An integer matrix dictating the layout of the plot;
#'   the default is a one-column matrix with one element for each item
#' @param main_title A character vector giving the plots' main titles;
#'   default is "Item Response Function".
#' @param sub An optional character vector of subtitles for the resulting plots,
#'   to be pasted onto the main title (helpful for titling individual plots
#'   when plotting multiple items' IRFs).
#' @param option_names An optional character vector giving names for the items'
#'   options; if NULL, generic names (e.g. "Option 1", "Option 2", etc.)
#'   are used
#' @param line_types An optional integer vector specifying lty for each option;
#'   if not provided, the first option for each question will have lty = 1,
#'   the second will have lty = 2, etc.
#' @param color A specification of the colors to draw lines in. Colors can be
#'   specified by either a character vector of colors (either names that R
#'   recognizes or hexidecimal specifications) or a function taking a single
#'   argument for the number of colors to return. See
#'   \code{\link{color_palettes}} for a list of color palettes provided
#'   by \code{bggum}. The default is "black" (for all lines).
#' @param rug A logical vector of length one specifying whether to draw a
#'   rug of theta estimates; the default is FALSE
#' @param thetas An optional vector of theta estimates for rug drawing;
#'   if \code{rug = TRUE} and \code{thetas} is not provided,
#'   an error will be thrown.
#' @param responses An optional matrix or vector (if the IRF for only one item
#'   is desired) of responses; if \code{rug = TRUE} and \code{responses} is not
#'   provided, an error will be thrown.
#'   NOTE: The lowest response for each item should be 0, not 1.
#' @param sides A vector giving the side(s) to draw the rug on if
#'   \code{rug = TRUE}; if the vector is of length > 1, the first option for
#'   each item will be drawn on the side given by the first element of the
#'   vector, the rug for the second option for each item will be drawn on the
#'   side given by the second element of the vector, etc.
#' @param rug_colors A vector giving the color(s) to draw the rug in if
#'   \code{rug = TRUE} or a function taking a single argument for the number of
#'   colors to return. See \code{\link{color_palettes}} for a list of color
#'   palettes provided by \code{bggum}. The default is "black" (for all rugs).
#'
#' @export
irf <- function(a, d, t, from = -3, to = 3, by = 0.01, layout_matrix = 1,
                main_title = "Item Response Function", sub = "",
                option_names = NULL, line_types = NULL, color = "black",
                rug = FALSE, thetas = NULL, responses = NULL, sides = 1,
                rug_colors = "black") {
    if ( length(a) != length(d) | (length(a) > 1 & length(a) != length(t)) ) {
        stop("Please provide a, d, and t of the same length.", call. = FALSE)
    }
    m <- length(a)
    if ( rug & is.null(responses) ) {
        stop("Please provide responses when rug = TRUE.", call. = FALSE)
    }
    if ( rug & is.null(thetas) ) {
        stop("Please provide theta estimates when rug = TRUE.", call. = FALSE)
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
            stop("Please provide responses for all items when rug = TRUE.",
                 call. = FALSE)
        }
    }
    if ( length(main_title) == 1 ) {
        main_title <- rep(main_title, m)
    }
    if ( length(sub) == 1 ) {
        sub <- rep(sub, m)
    }
    margins <- "if"(any(main_title != ""), c(3, 3, 3, 1), c(3, 3, 1, 1)) + 0.1
    opar <- graphics::par(mar = margins)
    on.exit(graphics::par(opar))
    graphics::layout(layout_matrix)
    th <- seq(from = from, to = to, by = by)
    if ( is.list(t) ) {
        K <- sapply(t, length)
    }
    else {
        K <- length(t)
        t <- list(t)
    }
    if ( is.character(color) ) {
        colors <- rep(color, length.out = max(K))
    } else if ( is.function(color) ) {
        colors <- color(max(K))
    } else {
        stop("Invalid color specification.")
    }
    option_vector <- 1:max(K)
    if ( is.null(option_names) ) {
        option_names <- paste("Option", option_vector)
    }
    if ( is.null(line_types) ) {
        line_types <- option_vector
    }
    else {
        if ( length(line_types) == 1 ) {
            line_types = rep(line_types, max(K))
        }
    }
    if ( length(sides) == 1 ) {
        sides <- rep(sides, max(K))
    }
    if ( is.character(rug_colors) ) {
        rug_colors <- rep(rug_colors, length.out = max(K))
    } else if ( is.function(rug_colors) ) {
        rug_colors <- rug_colors(max(K))
    } else {
        stop("Invalid rug color specification.")
    }
    for ( j in 1:m ) {
        y <- sapply(th, function(x) ggumProbability(0, x, a[j], d[j], t[[j]]))
        graphics::plot(x = th, y = y, col = colors[1], type = "l",
                       main = paste(main_title[j], sub[j]),
                       xlab = "", ylab = "", yaxt = "n", xaxt = "n",
                       xlim = c(from, to), ylim=c(0, 1.2),
                       lty = line_types[1])
        graphics::axis(side = 1, tick = FALSE, line = -0.75)
        graphics::axis(side = 2, tick = FALSE, line = -0.75,
                       at = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1"))
        graphics::title(ylab = expression(P[ij](k)),
                        xlab = expression(theta),
                        line = 1.5)
        for ( k in 2:K[j] ) {
            graphics::lines(th, sapply(th, function(x) {
                ggumProbability(k-1, x, a[j], d[j], t[[j]])
            }), lty = line_types[k], col = colors[k])
        }
        graphics::legend(x = from, y = 1.2, option_names, lty = line_types,
                         horiz = TRUE, bty = "n", col = colors)
        if ( rug ) {
            for ( k in 1:K[j] ) {
                idx <- which(responses[ , j] == k - 1)
                k_thetas <- thetas[idx]
                graphics::rug(k_thetas, side = sides[k], col = rug_colors[k])
            }
        }
    }
    graphics::layout(1) # resets layout
}
