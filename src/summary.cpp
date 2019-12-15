#include <Rcpp.h>

using namespace Rcpp;

// Note this is quantile type 8, not the R default of 7
// (type 8 is recommeded in Hyndman and Fan 1996)
// Also note that x is assumed to be sorted
double quantile(NumericVector x, double p) {
    int n = x.size();
    double h = ((n + 1.0/3.0) * p) - 2.0/3.0;
    int lo = (int)h;
    int hi = lo + 1;
    double lambda = h - lo;
    return x[lo] + lambda * (x[hi] - x[lo]);
}

NumericVector summarize_vector(NumericVector x) {
    NumericVector y = x.sort();
    NumericVector result(5);
    result[0] = quantile(y, 0.025);
    result[1] = quantile(y, 0.5);
    result[2] = mean(x);
    result[3] = quantile(y, 0.975);
    result[4] = sd(x);
    return result;
}

// [[Rcpp::export]]
NumericMatrix summarize_matrix(NumericMatrix x) {
    int m = x.ncol();
    NumericMatrix result(m, 5);
    for ( int i = 0; i < m; ++i ) {
        result(i, _) = summarize_vector(x(_, i));
    }
    return result;
}

