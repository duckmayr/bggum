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
NumericVector probRow(IntegerVector responseVec, double theta,
                      NumericVector alphas, NumericVector deltas, List taus);
NumericVector probCol(IntegerVector responseVec, NumericVector thetas,
                      double alpha, double delta, NumericVector taus);
NumericVector ggumProbability(NumericVector k, double theta, double alpha,
        double delta, NumericVector tau);
// pv stands for proposed value, cv for current value
double acceptanceTheta(IntegerVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD);
double acceptanceThetaNeg(IntegerVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD);
double acceptanceThetaPos(IntegerVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD);
double acceptanceAlpha(IntegerVector responses, NumericVector thetas,
        double cv, double delta, NumericVector taus, double SD);
double acceptanceDelta(IntegerVector responses, NumericVector thetas,
        double alpha, double cv, NumericVector taus, double SD);
double acceptanceTau(int k, IntegerVector responses, NumericVector thetas,
        double alpha, double delta, NumericVector taus, double SD);
NumericMatrix ggumMCMC(IntegerMatrix responseMatrix, IntegerVector Kvector,
                       int iterations, int low, int high);
double updateTheta(double cv, IntegerVector choices, NumericVector a,
        NumericVector d, List t, double temp);
double updateAlpha(double cv, IntegerVector choices, NumericVector th,
        double d, NumericVector t, double temp);
double updateDelta(double cv, IntegerVector choices, NumericVector th,
        double a, NumericVector t, double temp);
double updateTau(int k, IntegerVector choices, NumericVector th,
        double a, double d, NumericVector t, double temp);


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
double r_lst(double df, double mu, double sigma);


# endif
