#include "../inst/include/ggum.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector init_thetas(int n, double mean, double sd) {
    return rnorm(n, mean, sd);
}

// [[Rcpp::export]]
NumericVector init_alphas(int m, double shape1, double shape2, double a,
        double b) {
    return r4beta(m, shape1, shape2, a, b);
}

// [[Rcpp::export]]
NumericVector init_deltas(int m, double shape1, double shape2, double a,
        double b) {
    return r4beta(m, shape1, shape2, a, b);
}

// [[Rcpp::export]]
List init_taus(int m, double shape1, double shape2, double a, double b,
        IntegerVector K) {
    List taus(m);
    for ( int j = 0; j < m; ++j ){
        NumericVector thisTau(K[j]);
        thisTau[Range(1, K[j]-1)] = r4beta(K[j]-1, 2, 2, -2, 0);
        taus[j] = thisTau;
    }
    return taus;
}
