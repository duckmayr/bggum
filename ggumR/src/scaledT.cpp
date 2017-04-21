#include <Rcpp.h>

using namespace Rcpp;

//' The location-scale Student's T distribution.
//'
//' Density, distribution function, quantile function, and random number
//' generation for the Student's T distribution shifted by location parameter
//' 'mu' and scaled by 'sigma'.
//'
//' @param n A numeric vector of length one; the number of numbers to generate
//' @param x,q A numeric vector of the quantiles of interest
//' @param p A numeric vector of the probabilities of interest
//' @param df A numeric vector of length one; the degrees of freedom of the
//'   distribution
//' @param sigma A numeric vector of length one; the scale parameter
//' @param mu A numeric vector of length one; the shifting parameter
//'
//' @return 'dScaledT' gives the density, 'pScaledT' gives the distribution
//'     function, 'qScaledT' gives the quantile function, and 'rScaledT'
//'     generates random numbers from the distribution
//'     
//' @references Jackman, Simon. 2009. \emph{Bayesian Analysis for the Social
//'     Sciences}. Wiley. p. 507.
//'
//' @name Location-Scale T Distribution
//' @rdname scaledT-dist
void scaledT_distribution();

//' @rdname scaledT-dist
//' @export
//[[Rcpp::export]]
NumericVector dScaledT(NumericVector x, double df, double mu, double sigma){
    return (1 / sigma) * dt((x - mu)/sigma, df);
}

//' @rdname scaledT-dist
//' @export
//[[Rcpp::export]]
NumericVector pScaledT(NumericVector q, double df, double mu, double sigma){
    return pt((q - mu)/sigma, df);
}

//' @rdname scaledT-dist
//' @export
//[[Rcpp::export]]
NumericVector qScaledT(NumericVector p, double df, double mu, double sigma){
    return qt(p, df) * sigma + mu;
}

//' @rdname scaledT-dist
//' @export
//[[Rcpp::export]]
NumericVector rScaledT(int n, double df, double mu, double sigma){
    return rt(n, df) * sigma + mu;
}

//' @rdname scaledT-dist
//' @export
//[[Rcpp::export]]
double r_TruncNorm(double mu, double sigma, double a, double b){
    double F_a = R::pnorm(a, mu, sigma, 1, 0);
    double F_b = R::pnorm(b, mu, sigma, 1, 0);
    double q = R::qnorm(F_a + R::unif_rand() * (F_b - F_a), mu, sigma, 1, 0);
    return std::min(std::max(a, q), b);
}
