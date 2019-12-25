#include <Rcpp.h>

using namespace Rcpp;

//[[Rcpp::export]]
double prob(const int choice, const double th, const double a,
        const double d, const NumericVector& t){
    int K = t.size();
    double result = 0, numerator = 0, denominator = 0, tSum = 0;
    for ( int k = 0; k < K; ++k ) {
        tSum += t[k];
        numerator = exp(a * (k * (th - d) - tSum));
        numerator += exp(a * ((2*K - 1 - k) * (th - d) - tSum));
        denominator += numerator;
        if ( k == choice ) {
            result = numerator;
        }
    }
    return result / denominator;
}

//[[Rcpp::export]]
NumericVector probCol(const IntegerVector& choices, const NumericVector& thetas,
        const double a, const double d, const NumericVector& t){
    int n = choices.size();
    int K = t.size();
    NumericVector result(n);
    for ( int i = 0; i < n; ++i ){
        if ( IntegerVector::is_na(choices[i]) ) {
            result[i] = 1.0;
            continue;
        }
        double numerator = 0, denominator = 0, tSum = 0;
        for ( int k = 0; k < K; ++k ) {
            tSum += t[k];
            numerator = exp(a * (k * (thetas[i] - d) - tSum));
            numerator += exp(a * ((2*K - 1 - k) * (thetas[i] - d) - tSum));
            denominator += numerator;
            if ( k == choices[i] ) {
                result[i] = numerator;
            }
        }
        result[i] /= denominator;
    }
    return result;
}

//[[Rcpp::export]]
NumericVector probRow(const IntegerVector& choices, const double th,
        const NumericVector& a, const NumericVector& d, const List& t){
    int m = choices.size();
    NumericVector result(m);
    for ( int j = 0; j < m; ++j ){
        if ( IntegerVector::is_na(choices[j]) ) {
            result[j] = 1.0;
            continue;
        }
        double numerator = 0, denominator = 0, tSum = 0;
        double a_j = a[j], d_j = d[j];
        NumericVector t_j = as<NumericVector>(t[j]);
        int K = t_j.size();
        for ( int k = 0; k < K; ++k ) {
            tSum += t_j[k];
            numerator = exp(a_j * (k * (th - d_j) - tSum));
            numerator += exp(a_j * ((2*K - 1 - k) * (th - d_j) - tSum));
            denominator += numerator;
            if ( k == choices[j] ) {
                result[j] = numerator;
            }
        }
        result[j] /= denominator;
    }
    return result;
}

