# bggum

[![Travis-CI Build Status](https://travis-ci.org/duckmayr/bggum.svg?branch=master)](https://travis-ci.org/duckmayr/bggum)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/duckmayr/bggum?branch=master&svg=true)](https://ci.appveyor.com/project/duckmayr/bggum)
[![Coverage status](https://codecov.io/gh/duckmayr/bggum/branch/master/graph/badge.svg)](https://codecov.io/github/duckmayr/bggum?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/bggum)](https://cran.r-project.org/package=bggum)
[![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/last-month/bggum)](https://cranlogs.r-pkg.org/badges/last-month/bggum)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)


**bggum** provides an implementation for the Generalized Graded Unfolding Model (GGUM) of Roberts, Donoghue, and Laughlin (2000) for R.
Rather than provide tools for marginal maximum likelihood estimation, which is available from the package `GGUM` currently available [on CRAN](https://CRAN.R-project.org/package=GGUM), we provide an implementation of the MCMC algorithm in de la Torre et al (2006) as well as a Metropolis-coupled MCMC algorithm, both coded in C++ to allow for reasonable execution time using [Rcpp](https://github.com/RcppCore/Rcpp) (and using probability distribution functions written in C++ from [RcppDist](https://github.com/duckmayr/RcppDist)).


## Installation

You can install `bggum` from CRAN using

```r
install.packages("bggum")
```

Or, you can install the development version from GitHub via

```r
devtools::install_github("duckmayr/bggum")
```


## Usage

`bggum` provides tools for Bayesian estimation of GGUM parameters.
The package has a vignette with a reasonably in-depth practical guide to Bayesian estimation of GGUM parameters using `bggum`, so here we simply provide some highlights of available features:
 - Hyper-parameter tuning functions `tune_proposals()` to tune proposal densities (useful for either sampler) and `tune_temperatures()` to tune temperature schedules for the Metropolis-coupled MCMC (MC3) sampler
 - An MCMC (`ggumMCMC()`) and MC3 (`ggumMC3()`) sampler. Note: For the reasons discussed in [Duck-Mayr and Montgomery (2019)](http://jbduckmayr.com/papers/ggum.pdf), we prefer the MC3 sampling approach.
 - Like other IRT models, the GGUM exhibits rotational invariance, so we provide a function (`post_process()`) to post-process posterior samples to restrict attention to the reflective mode of interest.
 - The results of `ggumMCMC()` and `ggumMC3()` are `"ggum"` class objects, which have a `summary()` method giving parameter estimates and posterior quantile and standard deviation statistics.
 - We provide item response function (`irf()`) and item characteristic curve (`icc()`) plotting functions.

The vignette with a full example and step by step instructions can be accessed via 

```r
vignette("bggum")
```


## Contributing

Before contributing, please consult the contributing guidelines in CONTRIBUTING.md.


## License

GPL (>= 2)

