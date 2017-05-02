context("probability")

test_that("ggumProbability gets the correct probability", {
  expect_equal(ggumProbability(1, 3, 1, 2, c(-2, -1, 1, 2)), 0.46361,
               tolerance=0.00001)
  expect_equal(ggumProbability(1, 1, 1, 1, c(-2, -1, 1, 2)), 0.20603,
               tolerance=0.00001)
})
