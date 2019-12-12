context("Sanity checks in plotting functions")

test_that("irf() sanity checks work", {
    expect_error(irf(1, 2:3, list(0:-1)),
                 "Please provide a, d, and t of the same length.")
    expect_error(irf(1, 2, list(0:-1), responses = 1, rug = TRUE),
                 "Please provide theta estimates when rug = TRUE.")
    expect_error(irf(1, 2, list(0:-1), rug = TRUE),
                 "Please provide responses when rug = TRUE")
    expect_error(irf(1, 2, list(0:-1), responses = array(1, dim = c(1, 1, 1))),
                 "Please provide a vector or rectangular data structure")
    expect_error(irf(1, 2, list(0:-1), responses = matrix(0, ncol = 2)),
                 "Please provide responses for all items when rug = TRUE")
})

test_that("icc() sanity checks work", {
    expect_error(icc(1, 2:3, list(0:-1)),
                 "Please provide a, d, and t of the same length.")
    expect_error(icc(1, 2, list(0:-1), responses = 1, plot_responses = TRUE),
                 "Please provide theta estimates when plot_responses = TRUE.")
    expect_error(icc(1, 2, list(0:-1), plot_responses = TRUE),
                 "Please provide responses when plot_responses = TRUE")
    expect_error(icc(1, 2, list(0:-1), responses = array(1, dim = c(1, 1, 1))),
                 "Please provide a vector or rectangular data structure")
    expect_error(icc(1, 2, list(0:-1), responses = matrix(0, ncol = 2)),
                 "Please provide responses for all items when plot_responses = TRUE")
})
