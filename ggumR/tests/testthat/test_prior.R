context("Test if the function can get prior probability correctly")

test_that("Pr(alpha) is greater than or equal to 0",
          expect_gte(getPriorAlpha(runif(1, min=0.25, max=4)),  0))
test_that("Pr(alpha) is less than or equal to 1",
          expect_lte(getPriorAlpha(runif(1, min=0.25, max=4)),  1))

test_that("Pr(delta) is greater than or equal to 0",
          expect_gte(getPriorDelta(runif(1, min=-5, max=5)),  0))
test_that("Pr(delta) is less than or equal to 1",
          expect_lte(getPriorDelta(runif(1, min=-5, max=5)),  1))

test_that("Pr(tau) is greater than or equal to 0",
          expect_gte(getPriorTaus(runif(1, min=-6, max=6)), 0))
test_that("Pr(tau) is less than or equal to 1",
          expect_lte(getPriorTaus(runif(1, min=-6, max=6)), 1))

test_that("Pr(theta) is greater than or equal to 0",
          expect_gte(getPriorTheta(rnorm(1)), 0))
test_that("Pr(theta) is less than or equal to 1",
          expect_lte(getPriorTheta(rnorm(1)), 1))
