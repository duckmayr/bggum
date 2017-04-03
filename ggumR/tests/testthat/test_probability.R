context("probability")

test_that("Get the correct probability", {
  expect_that(probability(1, 2, 3, 1, 2, 3),
              equals(0.675972863168057))
  expect_that(probability(1, 4, 1, 1, 1, c(-2, -1, 1, 2)),
              equals(0.206031909190019))
})