#ifndef __UTILITIES__
#define __UTILITIES__

#include <Rcpp.h>

using namespace Rcpp;

// Throughout all the c++ files, we attempt to stick to the convention that
// n means the number of respondents, m means the number of questions, i refers
// to a specific respondent, j refers to a specific question, and k refers to a
// specific option for some question j.



// DEFINE FUNCTIONS


// The Generalized Graded Unfolding Model (GGUM)
double prob(int k, double theta, double alpha, double delta, NumericVector tau);
NumericVector probRow(NumericVector responseVec, double theta,
                      NumericVector alphas, NumericVector deltas, List taus);
NumericVector probCol(NumericVector responseVec, NumericVector thetas,
                      double alpha, double delta, NumericVector taus);
NumericVector ggumProbability(NumericVector k, double theta, double alpha,
        double delta, NumericVector tau);
// pv stands for proposed value, cv for current value
double acceptanceTheta(NumericVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD);
double acceptanceThetaNeg(NumericVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD);
double acceptanceThetaPos(NumericVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD);
double acceptanceAlpha(NumericVector responses, NumericVector thetas,
        double cv, double delta, NumericVector taus, double SD);
double acceptanceDelta(NumericVector responses, NumericVector thetas,
        double alpha, double cv, NumericVector taus, double SD);
double acceptanceTau(int k, NumericVector responses, NumericVector thetas,
        double alpha, double delta, NumericVector taus, double SD);
NumericMatrix ggumMCMC(NumericMatrix responseMatrix, IntegerVector Kvector,
                       int iterations, int low, int high);


// The Truncated Normal Distribution
NumericVector dtruncnorm(NumericVector x, double mean, double sd,
        double a, double b);
NumericVector ptruncnorm(NumericVector q, double mean, double sd,
        double a, double b);
NumericVector qtruncnorm(NumericVector p, double mean, double sd,
        double a, double b);
NumericVector rtruncnorm(int n, double mean, double sd, double a, double b);
// The following scalar version won't be available to the user
// It's just for faster calculation in the MCMC sampler
double r_truncnorm(double mean, double sd, double a, double b);
double d_truncnorm(double x, double mean, double SD, double a, double b);


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


// The Location-Scale T Distribution
NumericVector dlst(NumericVector x, double df, double mu, double sigma);
NumericVector plst(NumericVector q, double df, double mu, double sigma);
NumericVector qlst(NumericVector p, double df, double mu, double sigma);
NumericVector rlst(int n, double df, double mu, double sigma);
// The following scalar versions won't be available to the user
// They're just for use by other c++ functions
double d_lst(double x, double df, double mu, double sigma);
double p_lst(double q, double df, double mu, double sigma);


// The Truncated Location-Scale T Distribution
NumericVector dtrunclst(NumericVector x, double df, double mu, double sigma,
        double a, double b);
NumericVector ptrunclst(NumericVector q, double df, double mu, double sigma,
        double a, double b);
NumericVector qtrunclst(NumericVector p, double df, double mu, double sigma,
        double a, double b);
NumericVector rtrunclst(int n, double df, double mu, double sigma,
        double a, double b);
// The following scalar version won't be available to the user
// It's just for use by other c++ functions
double r_trunclst(double df, double mu, double sigma, double a, double b);
double d_trunclst(double x, double df, double mu, double sigma, double a,
        double b);


# endif
