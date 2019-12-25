context("Tuning functions")

n <- 100
m <- 10
K <- 2
set.seed(123)
ggum_sim   <- ggum_simulation(n, m, K)
sds <- tune_proposals(ggum_sim$response_matrix, 100)
temps <- tune_temperatures(ggum_sim$response_matrix, n_temps = 3,
                           temp_tune_iterations = 100, n_draws = 50,
                           sd_tune_iterations = 100)

test_that("tune_proposals() produces expected output", {
    expect_length(sds, 4)
    expect_length(sds[[1]], n)
    expect_length(sds[[2]], m)
    expect_length(sds[[3]], m)
    expect_length(sds[[4]], m)
})

test_that("tune_temperatures() produces expected output", {
    expect_equal(temps, c(1, 0.7253554, 0.6013790))
})
