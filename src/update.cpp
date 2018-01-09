#include "../inst/include/ggum.h"

using namespace Rcpp;

// Parameter updating functions for MCMC:

double update_theta_MCMC(const IntegerVector& responses, const double cv,
        const NumericVector& alphas, const NumericVector& deltas,
        const List& taus, const double SD){
    // pv stands for proposed value, cv for current value
    // First we propose a new value for the current theta;
    // we use a truncated location-scale T distribution:
    // double pv = r_lst(1, cv, SD);
    double pv = R::rnorm(cv, SD);
    // Then we get the prior probabilities that theta is equal to
    // the current value (cvPrior) and proposed value (pvPrior):
    double pvPrior = R::dnorm(pv, 0.0, 1.0, 0);
    double cvPrior = R::dnorm(cv, 0.0, 1.0, 0);
    // The (log) likelihood of the data given the current and proposed thetas:
    double cvL = sum(log(probRow(responses, cv, alphas, deltas, taus)));
    double pvL = sum(log(probRow(responses, pv, alphas, deltas, taus)));
    // Then we calculate the acceptance rate:
    double acceptRate = exp(pvL - cvL) * (pvPrior/cvPrior);
    // If the acceptance rate is greater than one, not a number,
    // or larger than a uniform deviate, we accept the proposal:
    if ( acceptRate > 1 || R::runif(0.0, 1.0) < acceptRate) {
        return pv;
    }
    // Otherwise we reject the proposal
    return cv;
}

//[[Rcpp::export]]
double update_theta_neg_MCMC(IntegerVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD){
    double pv = R::rnorm(cv, SD);
    double pvPrior = R::dnorm(pv, 0.0, 1.0, 0);
    double cvPrior = R::dnorm(cv, 0.0, 1.0, 0);
    double cvL = sum(log(probRow(responses, cv, alphas, deltas, taus)));
    double pvL = sum(log(probRow(responses, pv, alphas, deltas, taus)));
    double acceptRate = exp(pvL - cvL) * (pvPrior/cvPrior);
    if ( acceptRate > 1.0 || R::runif(0.0, 1.0) < acceptRate) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double update_theta_pos_MCMC(IntegerVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD){
    double pv = R::rnorm(cv, SD);
    double pvPrior = R::dnorm(pv, 0.0, 1.0, 0);
    double cvPrior = R::dnorm(cv, 0.0, 1.0, 0);
    double cvL = sum(log(probRow(responses, cv, alphas, deltas, taus)));
    double pvL = sum(log(probRow(responses, pv, alphas, deltas, taus)));
    double acceptRate = exp(pvL - cvL) * (pvPrior/cvPrior);
    if ( acceptRate > 1.0 || R::runif(0.0, 1.0) < acceptRate) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double update_alpha_MCMC(const IntegerVector& responses,
        const NumericVector& thetas, const double cv, const double delta,
        const NumericVector& taus, const double SD){
    double pv = R::rnorm(cv, SD);
    double pvPrior = d_4beta(pv, 1.5, 1.5, 0.25, 4.0);
    if ( pvPrior == 0.0 ) {
        return cv;
    }
    double cvPrior = d_4beta(cv, 1.5, 1.5, 0.25, 4.0);
    double cvL = sum(log(probCol(responses, thetas, cv, delta, taus)));
    double pvL = sum(log(probCol(responses, thetas, pv, delta, taus)));
    double acceptRate = exp(pvL - cvL) * (pvPrior/cvPrior);
    if ( acceptRate > 1.0 || R::runif(0.0, 1.0) < acceptRate) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double update_delta_MCMC(const IntegerVector& responses,
        const NumericVector& thetas, const double alpha, const double cv,
        const NumericVector& taus, const double SD){
    double pv = R::rnorm(cv, SD);
    double pvPrior = d_4beta(pv, 2.0, 2.0, -5.0, 5.0);
    if ( pvPrior == 0.0 ) {
        return cv;
    }
    double cvPrior = d_4beta(cv, 2.0, 2.0, -5.0, 5.0);
    double cvL = sum(log(probCol(responses, thetas, alpha, cv, taus)));
    double pvL = sum(log(probCol(responses, thetas, alpha, pv, taus)));
    double acceptRate = exp(pvL - cvL) * (pvPrior/cvPrior);
    if ( acceptRate > 1.0 || R::runif(0.0, 1.0) < acceptRate) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double update_tau_MCMC(const int k, const IntegerVector& responses,
        const NumericVector& thetas, const double alpha, const double delta,
        const NumericVector& taus, const double SD){
    // For taus, we need a copy of the entire vector
    NumericVector pv = clone(taus);
    // So that when we replace one of its values with a proposed value,
    // we can still compute probability of responses according to the GGUM
    // under the original tau vector
    pv[k] = R::rnorm(taus[k], SD);
    double pvPrior = d_4beta(pv[k], 2.0, 2.0, -6.0, 6.0);
    if ( pvPrior == 0.0 ) {
        return taus[k];
    }
    double cvPrior = d_4beta(taus[k], 2.0, 2.0, -6.0, 6.0);
    double cvL = sum(log(probCol(responses, thetas, alpha, delta, taus)));
    double pvL = sum(log(probCol(responses, thetas, alpha, delta, pv)));
    double acceptRate = exp(pvL - cvL) * (pvPrior/cvPrior);
    if ( acceptRate > 1.0 || R::runif(0.0, 1.0) < acceptRate) {
        return pv[k];
    }
    return taus[k];
}

// Parameter updating functions for MC3:

//[[Rcpp::export]]
double update_theta_MC3(const double cv, const IntegerVector& choices,
        const NumericVector& a, const NumericVector& d, const List& t,
        const double temp, const double SD){
    double pv = R::rnorm(cv, SD);
    double cvPrior = R::dnorm(cv, 0.0, 1.0, 0);
    double pvPrior = R::dnorm(pv, 0.0, 1.0, 0);
    double cvL = sum(log(probRow(choices, cv, a, d, t)));
    double pvL = sum(log(probRow(choices, pv, a, d, t)));
    double r = pow((exp(pvL - cvL)) * (pvPrior/cvPrior), temp);
    if ( r > 1.0 || R::runif(0.0, 1.0) < r) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double update_alpha_MC3(const double cv, const IntegerVector& choices,
        const NumericVector& th, const double d, const NumericVector& t,
        const double temp, const double SD){
    double pv = R::rnorm(cv, SD);
    double pvPrior = d_4beta(pv, 1.5, 1.5, 0.25, 4.0);
    if ( pvPrior == 0.0 ) {
        return cv;
    }
    double cvPrior = d_4beta(cv, 1.5, 1.5, 0.25, 4.0);
    double cvL = sum(log(probCol(choices, th, cv, d, t)));
    double pvL = sum(log(probCol(choices, th, pv, d, t)));
    double r = pow((exp(pvL - cvL)) * (pvPrior/cvPrior), temp);
    if ( r > 1.0 || R::runif(0.0, 1.0) < r) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double update_delta_MC3(const double cv, const IntegerVector& choices,
        const NumericVector& th, const double a, const NumericVector& t,
        const double temp, const double SD){
    double pv = R::rnorm(cv, SD);
    double pvPrior = d_4beta(pv, 2.0, 2.0, -5.0, 5.0);
    if ( pvPrior == 0.0 ) {
        return cv;
    }
    double cvPrior = d_4beta(cv, 2.0, 2.0, -5.0, 5.0);
    double cvL = sum(log(probCol(choices, th, a, cv, t)));
    double pvL = sum(log(probCol(choices, th, a, pv, t)));
    double r = pow((exp(pvL - cvL)) * (pvPrior/cvPrior), temp);
    if ( r > 1.0 || R::runif(0.0, 1.0) < r) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double update_tau_MC3(const int k, const IntegerVector& choices,
        const NumericVector& th, const double a, const double d,
        const NumericVector& t, const double temp, const double SD){
    NumericVector pv = clone(t);
    pv[k] = R::rnorm(t[k], SD);
    double pvPrior = d_4beta(pv[k], 2.0, 2.0, -6.0, 6.0);
    if ( pvPrior == 0.0 ) {
        return t[k];
    }
    double cvPrior = d_4beta(t[k], 2.0, 2.0, -6.0, 6.0);
    double cvL = sum(log(probCol(choices, th, a, d, t)));
    double pvL = sum(log(probCol(choices, th, a, d, pv)));
    double r = pow((exp(pvL - cvL)) * (pvPrior/cvPrior), temp);
    if ( r > 1.0 || R::runif(0.0, 1.0) < r) {
        return pv[k];
    }
    return t[k];
}

