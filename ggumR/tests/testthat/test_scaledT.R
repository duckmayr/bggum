context('Test if scaledT functions are accurate')

# Define some parameters for the test:
x <- seq(from=-1, to=5, by=0.05)
x2 <- seq(from=0, to=1, by=0.05)
df <- 1; mu <- 2; sigma <- 1

# Test c++ functions against manual R calculations:
test_that('dScaledT matches manual calculation',
          expect_equal(dScaledT(x, df, mu, sigma),
                       (1/sigma) * dt((x-mu)/sigma, df)))
test_that('pScaledT matches manual calculation',
          expect_equal(pScaledT(x, df, mu, sigma),
                       pt((x-mu)/sigma, df)))
test_that('qScaledT matches manual calculation',
          expect_equal(qScaledT(x2, df, mu, sigma),
                       (qt(x2, df) * sigma) + mu))

