context("Test if likelihoodCol is always between 0 and 1")
# Define the parameters to use in the test:
set.seed(123)
response_matrix <- cbind(sample(1:4, 1, replace = TRUE), 
                         sample(1:4, 1, replace = TRUE), 
                         sample(1:4, 1, replace = TRUE)) 
alphas <- runif(3, 0, 2)
deltas <- runif(3, 0, 1)
thetas <- rnorm(1)
taus <- list(c(0, runif(4, 0, 1)), c(0, runif(4, 0, 1)), c(0, runif(4, 0, 1)))

# Test likelihoodCol \in [0,1]
test_that("Output is equal or larger than 0", 
          expect_gte(likelihoodCol(thetas, response_matrix, alphas, deltas, taus), 0))
test_that("Output is equal or less than 1",
          expect_lte(likelihoodCol(thetas, response_matrix, alphas, deltas, taus), 1))
