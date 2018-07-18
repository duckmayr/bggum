#include "../inst/include/bggum.h"

using namespace Rcpp;

// Parameter updating functions for MCMC:

double update_theta_MCMC(const IntegerVector& responses, const double cv,
        const NumericVector& alphas, const NumericVector& deltas,
        const List& taus, const double proposal_sd,
        const double prior_mean, const double prior_sd){
    // pv stands for proposed value, cv for current value
    // First we propose a new value for the current theta;
    // we use a normal distribution centered on the current value:
    double pv = R::rnorm(cv, proposal_sd);
    // Then we get the (log) prior probabilities that theta is equal to
    // the current value (cvPrior) and proposed value (pvPrior):
    double pvPrior = R::dnorm(pv, prior_mean, prior_sd, 1);
    double cvPrior = R::dnorm(cv, prior_mean, prior_sd, 1);
    // The (log) likelihood of the data given the current and proposed thetas:
    double cvL = sum(log_probRow(responses, cv, alphas, deltas, taus));
    double pvL = sum(log_probRow(responses, pv, alphas, deltas, taus));
    // Then we calculate the acceptance rate:
    double acceptRate = pvL - cvL + pvPrior - cvPrior;
    // If the acceptance rate is greater than one, not a number,
    // or larger than a uniform deviate, we accept the proposal:
    if ( acceptRate > 0.0 || log(R::runif(0.0, 1.0)) < acceptRate) {
        return pv;
    }
    // Otherwise we reject the proposal
    return cv;
}

//[[Rcpp::export]]
double update_alpha_MCMC(const IntegerVector& responses,
        const NumericVector& thetas, const double cv, const double delta,
        const NumericVector& taus, const double proposal_sd,
        const double shape1, const double shape2,
        const double a, const double b){
    double pv = R::rnorm(cv, proposal_sd);
    double pvPrior = d_4beta(pv, shape1, shape2, a, b, 1);
    if ( pvPrior == R_NegInf ) {
        return cv;
    }
    double cvPrior = d_4beta(cv, shape1, shape2, a, b, 1);
    double cvL = sum(log_probCol(responses, thetas, cv, delta, taus));
    double pvL = sum(log_probCol(responses, thetas, pv, delta, taus));
    double acceptRate = pvL - cvL + pvPrior - cvPrior;
    if ( acceptRate > 0.0 || log(R::runif(0.0, 1.0)) < acceptRate) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double update_delta_MCMC(const IntegerVector& responses,
        const NumericVector& thetas, const double alpha, const double cv,
        const NumericVector& taus, const double proposal_sd,
        const double shape1, const double shape2,
        const double a, const double b){
    double pv = R::rnorm(cv, proposal_sd);
    double pvPrior = d_4beta(pv, shape1, shape2, a, b, 1);
    if ( pvPrior == R_NegInf ) {
        return cv;
    }
    double cvPrior = d_4beta(cv, shape1, shape2, a, b, 1);
    double cvL = sum(log_probCol(responses, thetas, alpha, cv, taus));
    double pvL = sum(log_probCol(responses, thetas, alpha, pv, taus));
    double acceptRate = pvL - cvL + pvPrior - cvPrior;
    if ( acceptRate > 0.0 || log(R::runif(0.0, 1.0)) < acceptRate) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double update_tau_MCMC(const int k, const IntegerVector& responses,
        const NumericVector& thetas, const double alpha, const double delta,
        const NumericVector& taus, const double proposal_sd,
        const double shape1, const double shape2,
        const double a, const double b){
    // For taus, we need a copy of the entire vector
    NumericVector pv = clone(taus);
    // So that when we replace one of its values with a proposed value,
    // we can still compute probability of responses according to the GGUM
    // under the original tau vector
    pv[k] = R::rnorm(taus[k], proposal_sd);
    double pvPrior = d_4beta(pv[k], shape1, shape2, a, b, 1);
    if ( pvPrior == R_NegInf ) {
        return taus[k];
    }
    double cvPrior = d_4beta(taus[k], shape1, shape2, a, b, 1);
    double cvL = sum(log_probCol(responses, thetas, alpha, delta, taus));
    double pvL = sum(log_probCol(responses, thetas, alpha, delta, pv));
    double acceptRate = pvL - cvL + pvPrior - cvPrior;
    if ( acceptRate > 0.0 || log(R::runif(0.0, 1.0)) < acceptRate) {
        return pv[k];
    }
    return taus[k];
}

// Parameter updating functions for MC3:

//[[Rcpp::export]]
double update_theta_MC3(const double cv, const IntegerVector& choices,
        const NumericVector& a, const NumericVector& d, const List& t,
        const double temp, const double proposal_sd,
        const double prior_mean, const double prior_sd){
    double pv = R::rnorm(cv, proposal_sd);
    double cvPrior = R::dnorm(cv, prior_mean, prior_sd, 1);
    double pvPrior = R::dnorm(pv, prior_mean, prior_sd, 1);
    double cvL = sum(log_probRow(choices, cv, a, d, t));
    double pvL = sum(log_probRow(choices, pv, a, d, t));
    double r = temp * (pvL - cvL + pvPrior - cvPrior);
    if ( r > 0.0 || log(R::runif(0.0, 1.0)) < r) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double update_alpha_MC3(const double cv, const IntegerVector& choices,
        const NumericVector& th, const double d, const NumericVector& t,
        const double temp, const double proposal_sd,
        const double shape1, const double shape2,
        const double a, const double b){
    double pv = R::rnorm(cv, proposal_sd);
    double pvPrior = d_4beta(pv, shape1, shape2, a, b, 1);
    if ( pvPrior == 0.0 ) {
        return cv;
    }
    double cvPrior = d_4beta(cv, shape1, shape2, a, b, 1);
    double cvL = sum(log_probCol(choices, th, cv, d, t));
    double pvL = sum(log_probCol(choices, th, pv, d, t));
    double r = temp * (pvL - cvL + pvPrior - cvPrior);
    if ( r > 0.0 || log(R::runif(0.0, 1.0)) < r) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double update_delta_MC3(const double cv, const IntegerVector& choices,
        const NumericVector& th, const double a, const NumericVector& t,
        const double temp, const double proposal_sd,
        const double shape1, const double shape2,
        const double low, const double high){
    double pv = R::rnorm(cv, proposal_sd);
    double pvPrior = d_4beta(pv, shape1, shape2, low, high, 1);
    if ( pvPrior == 0.0 ) {
        return cv;
    }
    double cvPrior = d_4beta(cv, shape1, shape2, low, high, 1);
    double cvL = sum(log_probCol(choices, th, a, cv, t));
    double pvL = sum(log_probCol(choices, th, a, pv, t));
    double r = temp * (pvL - cvL + pvPrior - cvPrior);
    if ( r > 0.0 || log(R::runif(0.0, 1.0)) < r) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double update_tau_MC3(const int k, const IntegerVector& choices,
        const NumericVector& th, const double a, const double d,
        const NumericVector& t, const double temp, const double proposal_sd,
        const double shape1, const double shape2,
        const double low, const double high){
    NumericVector pv = clone(t);
    pv[k] = R::rnorm(t[k], proposal_sd);
    double pvPrior = d_4beta(pv[k], shape1, shape2, low, high, 1);
    if ( pvPrior == 0.0 ) {
        return t[k];
    }
    double cvPrior = d_4beta(t[k], shape1, shape2, low, high, 1);
    double cvL = sum(log_probCol(choices, th, a, d, t));
    double pvL = sum(log_probCol(choices, th, a, d, pv));
    double r = temp * (pvL - cvL + pvPrior - cvPrior);
    if ( r > 0.0 || log(R::runif(0.0, 1.0)) < r) {
        return pv[k];
    }
    return t[k];
}

