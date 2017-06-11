#include "../inst/include/ggumR.h"

using namespace Rcpp;

//' The Truncated Normal Distribution
//'
//' Provides probability density, cumulative distribution, quantile, and random
//' number generating functions for the truncated normal distribution.
//'
//' @param x,q A numeric vector of quantiles
//' @param p A numeric vector of probabilities
//' @param n The number of observations to generate
//' @param mean The mean of the distribution
//' @param SD The standard deviation of the distribution
//' @param a The lowest possible value of the distribution
//' @param b The highest possible value of the distribution
//'
//' @return \code{dtruncnorm} gives gives the density, \code{ptruncnorm} gives
//'   the distribution function, \code{qtruncnorm} gives the quantile function,
//'   and \code{rtruncnorm} generates random numbers from the distribution.
//'
//' @name Truncated Normal
//' @rdname truncnorm
void truncnorm();

//' @rdname truncnorm
//' @export
//[[Rcpp::export]]
NumericVector dtruncnorm(NumericVector x, double mean, double SD,
        double a, double b){
    double scale = R::pnorm(b, mean, SD, 1, 0) - R::pnorm(a, mean, SD, 1, 0);
    NumericVector result = dnorm(x, mean, SD) / scale;
    for ( int l = 0; l < x.size(); l++ ) {
        if ( x[l] < a || x[l] > b ) {
            result[l] = 0;
        }
    }
    return result;
}

//' @rdname truncnorm
//' @export
//[[Rcpp::export]]
NumericVector ptruncnorm(NumericVector q, double mean, double SD,
        double a, double b){
    int N = q.size();
    NumericVector result(N);
    double F_a = R::pnorm(a, mean, SD, 1, 0);
    double F_b = R::pnorm(b, mean, SD, 1, 0);
    for ( int l = 0; l < N; l++ ) {
        result[l] = std::max(a, std::min(q[l], b));
    }
    return (pnorm(result, mean, SD) - F_a) / (F_b - F_a);
}

//' @rdname truncnorm
//' @export
//[[Rcpp::export]]
NumericVector qtruncnorm(NumericVector p, double mean, double SD,
        double a, double b){
    double F_a = R::pnorm(a, mean, SD, 1, 0);
    double F_b = R::pnorm(b, mean, SD, 1, 0);
    return pmin(b, pmax(a, qnorm(F_a + p * (F_b - F_a), mean, SD)));
}

//' @rdname truncnorm
//' @export
//[[Rcpp::export]]
NumericVector rtruncnorm(int n, double mean, double SD, double a, double b){
    return qtruncnorm(runif(n), mean, SD, a, b);
}

// The following scalar versions won't be available to the user
// They're just for faster calculation in the MCMC sampler
//[[Rcpp::export]]
double r_truncnorm(double mean, double SD, double a, double b){
    double F_a = R::pnorm(a, mean, SD, 1, 0);
    double F_b = R::pnorm(b, mean, SD, 1, 0);
    double q = R::qnorm(F_a + R::unif_rand() * (F_b - F_a), mean, SD, 1, 0);
    return std::min(std::max(a, q), b);
}

//[[Rcpp::export]]
double d_truncnorm(double x, double mean, double SD, double a, double b){
    if ( x < a || x > b ) {
        return 0;
    }
    double scale = R::pnorm(b, mean, SD, 1, 0) - R::pnorm(a, mean, SD, 1, 0);
    return R::dnorm(x, mean, SD, 0) / scale;
}
