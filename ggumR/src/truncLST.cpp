#include "../inst/include/ggumR.h"

using namespace Rcpp;

//' The Truncated Location-Scale T distribution.
//'
//' Provides probability density, cumulative distribution, quantile, and random
//' number generating functions for the truncated location-scale T distribution.
//'
//' @param x,q A numeric vector of quantiles
//' @param p A numeric vector of probabilities
//' @param n The number of observations to generate
//' @param df A numeric vector of length one; the degrees of freedom of the
//'   distribution
//' @param mu A numeric vector of length one; the shifting parameter
//' @param sigma A numeric vector of length one; the scale parameter
//' @param a The lowest possible value of the distribution
//' @param b The highest possible value of the distribution
//'
//' @return \code{dtrunclst} gives gives the density, \code{ptrunclst} gives
//'   the distribution function, \code{qtrunclst} gives the quantile function,
//'   and \code{rtrunclst} generates random numbers from the distribution.
//'
//' @name Truncated Location-Scale T
//' @rdname trunclst
void trunclst();

//' @rdname trunclst
//' @export
//[[Rcpp::export]]
NumericVector dtrunclst(NumericVector x, double df, double mu, double sigma,
        double a, double b){
    double scale = p_lst(b, df, mu, sigma) - p_lst(a, df, mu, sigma);
    NumericVector result = dlst(x, df, mu, sigma) / scale;
    for ( int l = 0; l < x.size(); l++ ) {
        if ( x[l] < a || x[l] > b ) {
            result[l] = 0;
        }
    }
    return result;
}

//' @rdname trunclst
//' @export
//[[Rcpp::export]]
NumericVector ptrunclst(NumericVector q, double df, double mu, double sigma,
        double a, double b){
    int N = q.size();
    NumericVector result(N);
    double F_a = p_lst(a, df, mu, sigma);
    double F_b = p_lst(b, df, mu, sigma);
    for ( int l = 0; l < N; l++ ) {
        result[l] = std::max(a, std::min(q[l], b));
    }
    return (plst(result, df, mu, sigma) - F_a) / (F_b - F_a);
}

//' @rdname trunclst
//' @export
//[[Rcpp::export]]
NumericVector qtrunclst(NumericVector p, double df, double mu, double sigma,
        double a, double b){
    double F_a = p_lst(a, df, mu, sigma);
    double F_b = p_lst(b, df, mu, sigma);
    return pmin(b, pmax(a, qlst(F_a + p * (F_b - F_a), df, mu, sigma)));
}

//' @rdname trunclst
//' @export
//[[Rcpp::export]]
NumericVector rtrunclst(int n, double df, double mu, double sigma,
        double a, double b){
    return qtrunclst(runif(n), df, mu, sigma, a, b);
}

// The following scalar versions won't be available to the user;
// they're just for faster calculation in the MCMC sampler
//[[Rcpp::export]]
double r_trunclst(double df, double mu, double sigma, double a, double b){
    double F_a = R::pt((a - mu)/sigma, df, 1, 0);
    double F_b = R::pt((b - mu)/sigma, df, 1, 0);
    double q = R::qt(F_a + R::unif_rand() * (F_b - F_a), df, 1, 0) * sigma + mu;
    return std::min(std::max(a, q), b);
}

//[[Rcpp::export]]
double d_trunclst(double x, double df, double mu, double sigma, double a,
        double b){
    if ( x < a || x > b ) {
        return 0;
    }
    double scale = p_lst(b, df, mu, sigma) - p_lst(a, df, mu, sigma);
    return d_lst(x, df, mu, sigma) / scale;
}
