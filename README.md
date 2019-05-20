# bggum

[![Travis-CI Build Status](https://travis-ci.org/duckmayr/bggum.svg?branch=master)](https://travis-ci.org/duckmayr/bggum)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/duckmayr/bggum?branch=master&svg=true)](https://ci.appveyor.com/project/duckmayr/bggum)
[![Coverage status](https://codecov.io/gh/duckmayr/bggum/branch/master/graph/badge.svg)](https://codecov.io/github/duckmayr/bggum?branch=master)

**bggum** provides an implementation for the Generalized Graded Unfolding Model (GGUM) of Roberts, Donoghue, and Laughlin (2000) for R. Rather than provide tools for marginal maximum likelihood estimation, which is available from the package `GGUM` currently available [on CRAN](https://CRAN.R-project.org/package=GGUM), we provide an implementation of the MCMC algorithm in de la Torre et al (2006) as well as a Metropolis-coupled MCMC algorithm, both coded in C++ to allow for reasonable execution time using [Rcpp](https://github.com/RcppCore/Rcpp) (and using probability distribution functions written in C++ from [RcppDist](https://github.com/duckmayr/RcppDist)).

## Installation

`bggum` is not yet on CRAN.
You can install the development version from GitHub via

```r
devtools::install_github("duckmayr/bggum")
```

## Contributing

Before contributing, please consult the contributing guidelines in CONTRIBUTING.md.

## License

GPL (>= 2)

