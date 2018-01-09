#include "../inst/include/ggum.h"

using namespace Rcpp;

//' The Four Parameter Beta Distribution
//'
//' Provides probability density, cumulative distribution, quantile, and random
//' number generating functions for the four parameter beta distribution.
//'
//' @param x,q A numeric vector of quantiles
//' @param p A numeric vector of probabilities
//' @param n The number of observations to generate
//' @param shape1 Shape 1
//' @param shape2 Shape 2
//' @param a The lowest possible value of the distribution
//' @param b The highest possible value of the distribution
//'
//' @return \code{d4beta} gives gives the density, \code{p4beta} gives the
//'   distribution function, \code{q4beta} gives the quantile function, and
//'   \code{r4beta} generates random numbers from the distribution.
//'
//' @name Four Parameter Beta
//' @rdname fourParamBeta
void fourParamBeta();

//' @rdname fourParamBeta
//' @export
//[[Rcpp::export]]
NumericVector d4beta(NumericVector x, double shape1, double shape2,
        double a, double b){
    NumericVector result = dbeta((x - a) / (b - a), shape1, shape2) / (b - a);
    for ( int l = 0; l < x.size(); l++ ) {
        if ( x[l] < a || x[l] > b ) {
            result[l] = 0;
        }
    }
    return result;
}

//' @rdname fourParamBeta
//' @export
//[[Rcpp::export]]
NumericVector p4beta(NumericVector q, double shape1, double shape2,
        double a, double b){
    return pbeta((q - a) / (b - a), shape1, shape2);
}

//' @rdname fourParamBeta
//' @export
//[[Rcpp::export]]
NumericVector q4beta(NumericVector p, double shape1, double shape2,
        double a, double b){
    return (b - a) * qbeta(p, shape1, shape2) + a;
}

//' @rdname fourParamBeta
//' @export
//[[Rcpp::export]]
NumericVector r4beta(int n, double shape1, double shape2, double a, double b){
    return (b - a) * rbeta(n, shape1, shape2) + a;
}

// The following scalar versions won't be available to the user
// They're just for faster calculation in the MCMC sampler
//[[Rcpp::export]]
double d_4beta(double x, double shape1, double shape2, double a, double b){
    if ( x < a || x > b ) {
        return 0.0;
    }
    return R::dbeta((x - a) / (b - a), shape1, shape2, 0) / (b - a);
}
//[[Rcpp::export]]
double r_4beta(double shape1, double shape2, double a, double b){
    return (b - a) * R::rbeta(shape1, shape2) + a;
}
