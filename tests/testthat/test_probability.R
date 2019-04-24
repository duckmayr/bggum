context("ggumProbability")

test_that("ggumProbability gives consistent exact answers", {
  expect_equal(ggumProbability(0, 3, 1, 2, c(-2, -1, 1, 2)), 0.46361,
               tolerance=0.00001)
  expect_equal(ggumProbability(0, 1, 1, 1, c(-2, -1, 1, 2)), 0.20603,
               tolerance=0.00001)
})

test_that("ggumProbability results match de la Torre et al. (2006) Fig. 1", {
    alpha <- 2
    delta <- 0
    tau <- c(0, -1, -0.7, -0.4)
    K <- 4
    # For k = 3 and theta = 0
    expect_lte(ggumProbability(3, 0, alpha, delta, tau),  0.65)
    expect_gte(ggumProbability(3, 0, alpha, delta, tau),  0.6)
    # For k = 2 and theta = 0.25
    expect_lte(ggumProbability(2, 0.25, alpha, delta, tau),  0.35)
    expect_gte(ggumProbability(2, 0.25, alpha, delta, tau),  0.3)
    # For k = 1 and theta = 1
    expect_lte(ggumProbability(1, 1, alpha, delta, tau),  0.40)
    expect_gte(ggumProbability(1, 1, alpha, delta, tau),  0.35)
    # For k = 1 and theta = -3
    expect_lte(ggumProbability(0, -3, alpha, delta, tau),  1)
    expect_gte(ggumProbability(0, -3, alpha, delta, tau),  0.96)
})

n <- 100
m <- 10
K <- 2
set.seed(123)
ggum_sim <- ggum_simulation(n, m, K)

test_that("ggumProbability handles matrices correctly", {
    expect_error(ggumProbability(ggum_sim$response_matrix, 1, 1, 1, 1),
                 "For a response matrix")
    probs <- ggumProbability(ggum_sim$response_matrix, ggum_sim$theta,
                             ggum_sim$alpha, ggum_sim$delta, ggum_sim$tau)
    expect_true(all(probs >= 0))
    expect_true(all(probs <= 1))
    expect_equal(dim(probs), dim(ggum_sim$response_matrix))
})

test_that("ggumProbability handles rows correctly", {
    expect_error(ggumProbability(1, 1:2, 1:2, 1:2, 1:2),
                 "For multiple items and respondents")
    expect_error(ggumProbability(1, 1:2, 1, 1, 1),
                 "Provide a response for each respondent.")
    probs <- ggumProbability(ggum_sim$response_matrix[1, ], ggum_sim$theta[1],
                             ggum_sim$alpha, ggum_sim$delta, ggum_sim$tau)
    expect_true(all(probs >= 0))
    expect_true(all(probs <= 1))
    expect_length(probs, m)
})

test_that("ggumProbability handles columns correctly", {
    expect_error(ggumProbability(1, 1, 1:2, 1:2, 1:2),
                 "Provide a response for each item.")
    probs <- ggumProbability(ggum_sim$response_matrix[ , 1], ggum_sim$theta,
                             ggum_sim$alpha[1], ggum_sim$delta[1],
                             ggum_sim$tau[[1]])
    expect_true(all(probs >= 0))
    expect_true(all(probs <= 1))
    expect_length(probs, n)
})
