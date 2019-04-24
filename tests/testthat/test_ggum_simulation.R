context("GGUM Simulation")

n <- 100
m <- 10
K <- 2
set.seed(123)
ggum_sim <- ggum_simulation(n, m, K)

test_that("ggum_simulation output is structured as expected", {
    expect_length(ggum_sim, 5)
    expect_named(ggum_sim, c("theta", "alpha", "delta", "tau", "response_matrix"))
    expect_setequal(sapply(ggum_sim, class), c("numeric", "list", "matrix"))
    expect_length(ggum_sim$theta, n)
    expect_length(ggum_sim$alpha, m)
    expect_length(ggum_sim$delta, m)
    expect_length(ggum_sim$tau, m)
    expect_equal(dim(ggum_sim$response_matrix), c(n, m))
})

test_that("ggum_simulation output values are as expected", {
    expect_setequal(c(ggum_sim$response_matrix), 0:1)
    expect_setequal(sapply(ggum_sim$tau, '[[', 1), 0)
    expect_true(all(sapply(ggum_sim$tau, '[[', 2) <=  0))
    expect_true(all(sapply(ggum_sim$tau, '[[', 2) >=  -2))
    expect_true(all(ggum_sim$alpha >=  0.25))
    expect_true(all(ggum_sim$alpha <=  4))
    expect_false(any(is.na(unlist(ggum_sim))))
})
