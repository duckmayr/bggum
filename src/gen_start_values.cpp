#include <Rcpp.h>
#include <4beta.h>

// [[Rcpp::export]]
Rcpp::NumericVector init_thetas(int n, double mean, double sd) {
    return Rcpp::rnorm(n, mean, sd);
}

// [[Rcpp::export]]
Rcpp::NumericVector init_alphas(int m, double shape1, double shape2, double a,
                                double b) {
    return r4beta(m, shape1, shape2, a, b);
}

// [[Rcpp::export]]
Rcpp::NumericVector init_deltas(int m, double shape1, double shape2, double a,
                                double b) {
    return r4beta(m, shape1, shape2, a, b);
}

// [[Rcpp::export]]
Rcpp::List init_taus(int m, double shape1, double shape2, double a, double b,
                     Rcpp::IntegerVector K) {
    Rcpp::List taus(m);
    for ( int j = 0; j < m; ++j ){
        Rcpp::NumericVector this_tau(K[j]);
        this_tau[Rcpp::Range(1, K[j]-1)] = r4beta(K[j]-1, shape1, shape2, a, b);
        taus[j] = this_tau;
    }
    return taus;
}

