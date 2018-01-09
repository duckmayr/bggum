context("Test if likelihoodCol is always between 0 and 1")
# Define the parameters to use in the test:
set.seed(123)
response_vector <- sample(1:4, 3, replace = TRUE)
alphas <- runif(1, 0, 2)
deltas <- runif(1, 0, 1)
thetas <- rnorm(3)
taus <- c(0, runif(3, 0, 1))

# Test likelihoodCol \in [0,1]
test_that("Output is equal or larger than 0", 
          expect_gte(exp(loglikelihoodCol(response_vector, thetas, alphas,
                         deltas, taus)), 0))
test_that("Output is equal or less than 1",
          expect_lte(exp(loglikelihoodCol(response_vector, thetas, alphas,
                         deltas, taus)), 1))
