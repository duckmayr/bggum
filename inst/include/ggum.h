#ifndef GGUM_H
#define GGUM_H

#include <Rcpp.h>

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
NumericVector ggumProbability(NumericVector k, double theta, double alpha,
        double delta, NumericVector tau);
NumericVector log_probCol(const IntegerVector& choices,
        const NumericVector& thetas, const double a, const double d,
        const NumericVector& t);
NumericVector log_probRow(const IntegerVector& choices, const double th,
        const NumericVector& a, const NumericVector& d, const List& t);

// Updating functions (pv stands for proposed value, cv for current value):
// For MCMC:
double update_theta_MCMC(const IntegerVector& responses, const double cv,
        const NumericVector& alphas, const NumericVector& deltas,
        const List& taus, const double SD);
double update_theta_neg_MCMC(IntegerVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD);
double update_theta_pos_MCMC(IntegerVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD);
double update_alpha_MCMC(const IntegerVector& responses,
        const NumericVector& thetas, const double cv, const double delta,
        const NumericVector& taus, const double SD);
double update_delta_MCMC(const IntegerVector& responses,
        const NumericVector& thetas, const double alpha, const double cv,
        const NumericVector& taus, const double SD);
double update_tau_MCMC(const int k, const IntegerVector& responses,
        const NumericVector& thetas, const double alpha, const double delta,
        const NumericVector& taus, const double SD);
// For MC3:
double update_theta_MC3(const double cv, const IntegerVector& choices,
        const NumericVector& a, const NumericVector& d, const List& t,
        const double temp, const double SD);
double update_alpha_MC3(const double cv, const IntegerVector& choices,
        const NumericVector& th, const double d, const NumericVector& t,
        const double temp, const double SD);
double update_delta_MC3(const double cv, const IntegerVector& choices,
        const NumericVector& th, const double a, const NumericVector& t,
        const double temp, const double SD);
double update_tau_MC3(const int k, const IntegerVector& choices,
        const NumericVector& th, const double a, const double d,
        const NumericVector& t, const double temp, const double SD);
// Proposal tuning function
List tune_proposals(const IntegerMatrix& responseMatrix, NumericVector& thetas,
        NumericVector& alphas, NumericVector& deltas, List& taus,
        const IntegerVector& K, const int burn_iters, int n, int m);

// The Four Parameter Beta Distribution
NumericVector d4beta(NumericVector x, double shape1, double shape2,
        double a, double b);
NumericVector p4beta(NumericVector q, double shape1, double shape2,
        double a, double b);
NumericVector q4beta(NumericVector p, double shape1, double shape2,
        double a, double b);
NumericVector r4beta(int n, double shape1, double shape2, double a, double b);
// The following scalar versions won't be available to the user
// They're just for faster calculation in the MCMC sampler
double d_4beta(double x, double shape1, double shape2, double a, double b);
double r_4beta(double shape1, double shape2, double a, double b);


# endif
