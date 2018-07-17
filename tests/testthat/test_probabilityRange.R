context("Test if probability produces results in a range of values that is
        similar to the one observed in de la Torres 2006")
# Define parameter to use in the test.
alpha <- 2
delta <- 0
tau <- c(0, -1, -0.7, -0.4)
K <- 4

# For k = 4 and theta = 0
test_that("Output is less than 0.65",
          expect_lte(ggumProbability(3, 0, alpha, delta, tau),  0.65))
test_that("Output is greater than 0.6",
          expect_gte(ggumProbability(3, 0, alpha, delta, tau),  0.6))

# For k = 3 and theta = 0.25
test_that("Output is less than 0.35",
          expect_lte(ggumProbability(2, 0.25, alpha, delta, tau),  0.35))
test_that("Output is greater than 0.3",
          expect_gte(ggumProbability(2, 0.25, alpha, delta, tau),  0.3))

# For k = 2 and theta = 1
test_that("Output is less than 0.40",
          expect_lte(ggumProbability(1, 1, alpha, delta, tau),  0.40))
test_that("Output is greater than 0.35",
          expect_gte(ggumProbability(1, 1, alpha, delta, tau),  0.35))

# For k = 1 and theta = -3
test_that("Output is less than 1",
          expect_lte(ggumProbability(0, -3, alpha, delta, tau),  1))
test_that("Output is greater than 0.96",
          expect_gte(ggumProbability(0, -3, alpha, delta, tau),  0.96))
