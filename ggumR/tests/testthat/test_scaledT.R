context('Test if location-scale T distribution functions are accurate')

# Define some parameters for the test:
x <- seq(from=-1, to=5, by=0.05)
x2 <- seq(from=0, to=1, by=0.05)
df <- 1; mu <- 2; sigma <- 1

# Test c++ functions against manual R calculations:
test_that('dlst matches manual calculation',
          expect_equal(dlst(x, df, mu, sigma),
                       (1/sigma) * dt((x-mu)/sigma, df)))
test_that('plst matches manual calculation',
          expect_equal(plst(x, df, mu, sigma),
                       pt((x-mu)/sigma, df)))
test_that('qlst matches manual calculation',
          expect_equal(qlst(x2, df, mu, sigma),
                       (qt(x2, df) * sigma) + mu))

