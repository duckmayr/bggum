#include "bggum.h"

using namespace Rcpp;

NumericVector sim_probs(const double th, const double a, const double d,
        const NumericVector& t){
    int K = t.size();
    NumericVector result(K);
    double numerator = 0, denominator = 0, tSum = 0;
    for ( int k = 0; k < K; ++k ) {
        tSum += t[k];
        numerator = exp(a * (k * (th - d) - tSum));
        numerator += exp(a * ((2*K - 1 - k) * (th - d) - tSum));
        denominator += numerator;
        result[k] = numerator;
    }
    return cumsum(result / denominator);
}

//[[Rcpp::export(.ggum_simulation)]]
IntegerMatrix ggum_simulation(const int n, const int m, const IntegerVector& K,
        const NumericVector& theta, const NumericVector& alpha,
        const NumericVector& delta, const List& tau) {
    IntegerMatrix resp_mat(n, m);
    for ( int i = 0; i < n; ++i ) {
        for ( int j = 0; j < m; ++j ) {
            int kk = K[j];
            double u = R::runif(0.0, 1.0);
            NumericVector probs = sim_probs(theta[i], alpha[j], delta[j],
                as<NumericVector>(tau[j]));
            for ( int k = 0; k < kk; ++k ) {
                if ( u < probs[k] ) {
                    resp_mat(i, j) = k;
                    break;
                }
            }
        }
    }
    return resp_mat;
}
