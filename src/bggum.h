#ifndef BGGUM_H
#define BGGUM_H

#include <Rcpp.h>
#include <4beta.h>

using namespace Rcpp;

// Throughout all the c++ files, we attempt to stick to the convention that
// n means the number of respondents, m means the number of questions, i refers
// to a specific respondent, j refers to a specific question, and k refers to a
// specific option for some question j.



// DEFINE FUNCTIONS


// The Generalized Graded Unfolding Model (GGUM)
// Probability functions:
double prob(const int choice, const double th, const double a,
        const double d, const NumericVector& t);
NumericVector probCol(const IntegerVector& choices, const NumericVector& thetas,
        const double a, const double d, const NumericVector& t);
NumericVector probRow(const IntegerVector& choices, const double th,
        const NumericVector& a, const NumericVector& d, const List& t);
NumericVector log_probCol(const IntegerVector& choices,
        const NumericVector& thetas, const double a, const double d,
        const NumericVector& t);
NumericVector log_probRow(const IntegerVector& choices, const double th,
        const NumericVector& a, const NumericVector& d, const List& t);

// Updating functions (pv stands for proposed value, cv for current value):
// For MCMC:
double update_theta_MCMC(const IntegerVector& responses, const double cv,
        const NumericVector& alphas, const NumericVector& deltas,
        const List& taus, const double proposal_sd,
        const double prior_mean, const double prior_sd);
double update_alpha_MCMC(const IntegerVector& responses,
        const NumericVector& thetas, const double cv, const double delta,
        const NumericVector& taus, const double proposal_sd,
        const double shape1, const double shape2,
        const double a, const double b);
double update_delta_MCMC(const IntegerVector& responses,
        const NumericVector& thetas, const double alpha, const double cv,
        const NumericVector& taus, const double proposal_sd,
        const double shape1, const double shape2,
        const double a, const double b);
double update_tau_MCMC(const int k, const IntegerVector& responses,
        const NumericVector& thetas, const double alpha, const double delta,
        const NumericVector& taus, const double proposal_sd,
        const double shape1, const double shape2,
        const double a, const double b);
// For MC3:
double update_theta_MC3(const double cv, const IntegerVector& choices,
        const NumericVector& a, const NumericVector& d, const List& t,
        const double temp, const double proposal_sd,
        const double prior_mean, const double prior_sd);
double update_alpha_MC3(const double cv, const IntegerVector& choices,
        const NumericVector& th, const double d, const NumericVector& t,
        const double temp, const double proposal_sd,
        const double shape1, const double shape2,
        const double a, const double b);
double update_delta_MC3(const double cv, const IntegerVector& choices,
        const NumericVector& th, const double a, const NumericVector& t,
        const double temp, const double proposal_sd,
        const double shape1, const double shape2,
        const double low, const double high);
double update_tau_MC3(const int k, const IntegerVector& choices,
        const NumericVector& th, const double a, const double d,
        const NumericVector& t, const double temp, const double proposal_sd,
        const double shape1, const double shape2,
        const double low, const double high);
// Proposal tuning function
List tune_proposals(const IntegerMatrix& responseMatrix, NumericVector& thetas,
        NumericVector& alphas, NumericVector& deltas, List& taus,
        const IntegerVector& K, const int tune_iters, int n, int m,
        double th_prior_mean, double th_prior_sd, double a_shape1,
        double a_shape2, double a_a, double a_b, double d_shape1,
        double d_shape2, double d_a, double d_b, double t_shape1,
        double t_shape2, double t_a, double t_b);
        

# endif
