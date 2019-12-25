#include "model.h"

double Model::theta_ll(double val, int t, int i) {
    double result = 0.0;
    int M = respondent_non_na_indices[i].size();
    for ( int j = 0; j < M; ++j ){
        int idx = respondent_non_na_indices[i][j];
        double numerator = 0, denominator = 0, chosen_numerator = 0;
        double dist = val - delta(idx, t);
        int Kj = K[idx];
        Rcpp::NumericVector tau_cumsum = Rcpp::cumsum(tau[t][idx]);
        for ( int k = 0; k < Kj; ++k ) {
            double tau_sum = tau_cumsum[k];
            numerator = exp(alpha(idx, t) * (k * dist - tau_sum));
            numerator += exp(alpha(idx, t) * ((2*Kj - 1 - k) * dist - tau_sum));
            denominator += numerator;
            if ( k == data(i, idx) ) {
                chosen_numerator = numerator;
            }
        }
        result += (log(chosen_numerator) - log(denominator));
    }
    return result;
}

double Model::alpha_ll(double val, int t, int j) {
    double result = 0.0;
    int N = item_non_na_indices[j].size();
    Rcpp::NumericVector tau_cumsum = Rcpp::cumsum(tau[t][j]);
    for ( int i = 0; i < N; ++i ){
        int idx = item_non_na_indices[j][i];
        double numerator = 0, denominator = 0, chosen_numerator = 0;
        double dist = theta(idx, t) - delta(j, t);
        int Kj = K[j];
        for ( int k = 0; k < Kj; ++k ) {
            double tau_sum = tau_cumsum[k];
            numerator = exp(val * (k * dist - tau_sum));
            numerator += exp(val * ((2*Kj - 1 - k) * dist - tau_sum));
            denominator += numerator;
            if ( k == data(idx, j) ) {
                chosen_numerator = numerator;
            }
        }
        result += (log(chosen_numerator) - log(denominator));
    }
    return result;
}

double Model::delta_ll(double val, int t, int j) {
    double result = 0.0;
    int N = item_non_na_indices[j].size();
    Rcpp::NumericVector tau_cumsum = Rcpp::cumsum(tau[t][j]);
    for ( int i = 0; i < N; ++i ){
        double numerator = 0, denominator = 0, chosen_numerator = 0;
        int idx = item_non_na_indices[j][i];
        double dist = theta(idx, t) - val;
        int Kj = K[j];
        for ( int k = 0; k < Kj; ++k ) {
            double tau_sum = tau_cumsum[k];
            numerator = exp(alpha(j, t) * (k * dist - tau_sum));
            numerator += exp(alpha(j, t) * ((2*Kj - 1 - k) * dist - tau_sum));
            denominator += numerator;
            if ( k == data(idx, j) ) {
                chosen_numerator = numerator;
            }
        }
        result += (log(chosen_numerator) - log(denominator));
    }
    return result;
}

double Model::tau_ll(Rcpp::NumericVector val, int t, int j) {
    double result = 0.0;
    Rcpp::NumericVector tau_cumsum = Rcpp::cumsum(val);
    int N = item_non_na_indices[j].size();
    for ( int i = 0; i < N; ++i ){
        int idx = item_non_na_indices[j][i];
        double numerator = 0, denominator = 0, chosen_numerator = 0;
        double dist = theta(idx, t) - delta(j, t);
        int Kj = K[j];
        for ( int k = 0; k < Kj; ++k ) {
            double tau_sum = tau_cumsum[k];
            numerator = exp(alpha(j, t) * (k * dist - tau_sum));
            numerator += exp(alpha(j, t) * ((2*Kj - 1 - k) * dist - tau_sum));
            denominator += numerator;
            if ( k == data(idx, j) ) {
                chosen_numerator = numerator;
            }
        }
        result += (log(chosen_numerator) - log(denominator));
    }
    return result;
}

double Model::item_ll(int t, int j) {
    double result = 0.0;
    Rcpp::NumericVector tau_cumsum = Rcpp::cumsum(tau[t][j]);
    int N = item_non_na_indices[j].size();
    for ( int i = 0; i < N; ++i ){
        int idx = item_non_na_indices[j][i];
        double numerator = 0, denominator = 0, chosen_numerator = 0;
        double dist = theta(idx, t) - delta(j, t);
        int Kj = K[j];
        for ( int k = 0; k < Kj; ++k ) {
            double tau_sum = tau_cumsum[k];
            numerator = exp(alpha(j, t) * (k * dist - tau_sum));
            numerator += exp(alpha(j, t) * ((2*Kj - 1 - k) * dist - tau_sum));
            denominator += numerator;
            if ( k == data(idx, j) ) {
                chosen_numerator = numerator;
            }
        }
        result += (log(chosen_numerator) - log(denominator));
    }
    return result;
}

double Model::log_likelihood(int t) {
    double result = 0.0;
    for ( int j = 0; j < m; ++j ) {
        result += item_ll(t, j);
    }
    return result;
}

