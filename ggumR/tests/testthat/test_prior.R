context("Test if the function can get prior probability correctly")

test_that("alpha is greater than or equal to 0",
          expect_gte(getPrior(param = "alpha", value = c(1, 2, 3, 4)),  0))
test_that("alpha is less than or equal to 1",
          expect_lte(getPrior(param = "alpha", value = c(1, 2, 3, 4)),  1))

test_that("delta is greater than or equal to 0",
          expect_gte(getPrior(param = "delta", value = c(1, 2, 3, 4)),  0))
test_that("delta is less than or equal to 1",
          expect_lte(getPrior(param = "delta", value = c(1, 2, 3, 4)),  1))

test_that("tau is greater than or equal to 0",
          expect_gte(getPrior(param = "tau", value = c(1, 2, 3, 4)),  0))
test_that("tau is less than or equal to 1",
          expect_lte(getPrior(param = "tau", value = c(1, 2, 3, 4)),  1))

test_that("theta is greater than or equal to 0",
          expect_gte(getPrior(param = "theta", value = c(1, 2, 3, 4)),  0))
test_that("theta is less than or equal to 1",
          expect_lte(getPrior(param = "theta", value = c(1, 2, 3, 4)),  1))
