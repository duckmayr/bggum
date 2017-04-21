#include "../inst/include/ggumR.h"

using namespace Rcpp;

//' The location-scale T distribution.
//'
//' Density, distribution function, quantile function, and random number
//' generation for the T distribution shifted by location parameter 'mu' and
//' scaled by 'sigma'.
//'
//' @param n A numeric vector of length one; the number of numbers to generate
//' @param x,q A numeric vector of the quantiles of interest
//' @param p A numeric vector of the probabilities of interest
//' @param df A numeric vector of length one; the degrees of freedom of the
//'   distribution
//' @param sigma A numeric vector of length one; the scale parameter
//' @param mu A numeric vector of length one; the shifting parameter
//'
//' @return 'dlst' gives the density, 'plst' gives the distribution function,
//'     'qlst' gives the quantile function, and 'rlst' generates random numbers
//'     from the distribution
//'     
//' @references Jackman, Simon. 2009. \emph{Bayesian Analysis for the Social
//'     Sciences}. Wiley. p. 507.
//'
//' @name Location-Scale T
//' @rdname scaledT-dist
void scaledT_distribution();

//' @rdname scaledT-dist
//' @export
//[[Rcpp::export]]
NumericVector dlst(NumericVector x, double df, double mu, double sigma){
    return (1 / sigma) * dt((x - mu)/sigma, df);
}

//' @rdname scaledT-dist
//' @export
//[[Rcpp::export]]
NumericVector plst(NumericVector q, double df, double mu, double sigma){
    return pt((q - mu)/sigma, df);
}

//' @rdname scaledT-dist
//' @export
//[[Rcpp::export]]
NumericVector qlst(NumericVector p, double df, double mu, double sigma){
    return qt(p, df) * sigma + mu;
}

//' @rdname scaledT-dist
//' @export
//[[Rcpp::export]]
NumericVector rlst(int n, double df, double mu, double sigma){
    return rt(n, df) * sigma + mu;
}

