// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/ggum.h"
#include <Rcpp.h>

using namespace Rcpp;

// d4beta
NumericVector d4beta(NumericVector x, double shape1, double shape2, double a, double b);
RcppExport SEXP _ggum_d4beta(SEXP xSEXP, SEXP shape1SEXP, SEXP shape2SEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< double >::type shape2(shape2SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(d4beta(x, shape1, shape2, a, b));
    return rcpp_result_gen;
END_RCPP
}
// p4beta
NumericVector p4beta(NumericVector q, double shape1, double shape2, double a, double b);
RcppExport SEXP _ggum_p4beta(SEXP qSEXP, SEXP shape1SEXP, SEXP shape2SEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< double >::type shape2(shape2SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(p4beta(q, shape1, shape2, a, b));
    return rcpp_result_gen;
END_RCPP
}
// q4beta
NumericVector q4beta(NumericVector p, double shape1, double shape2, double a, double b);
RcppExport SEXP _ggum_q4beta(SEXP pSEXP, SEXP shape1SEXP, SEXP shape2SEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< double >::type shape2(shape2SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(q4beta(p, shape1, shape2, a, b));
    return rcpp_result_gen;
END_RCPP
}
// r4beta
NumericVector r4beta(int n, double shape1, double shape2, double a, double b);
RcppExport SEXP _ggum_r4beta(SEXP nSEXP, SEXP shape1SEXP, SEXP shape2SEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< double >::type shape2(shape2SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(r4beta(n, shape1, shape2, a, b));
    return rcpp_result_gen;
END_RCPP
}
// d_4beta
double d_4beta(double x, double shape1, double shape2, double a, double b);
RcppExport SEXP _ggum_d_4beta(SEXP xSEXP, SEXP shape1SEXP, SEXP shape2SEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< double >::type shape2(shape2SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(d_4beta(x, shape1, shape2, a, b));
    return rcpp_result_gen;
END_RCPP
}
// r_4beta
double r_4beta(double shape1, double shape2, double a, double b);
RcppExport SEXP _ggum_r_4beta(SEXP shape1SEXP, SEXP shape2SEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< double >::type shape2(shape2SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(r_4beta(shape1, shape2, a, b));
    return rcpp_result_gen;
END_RCPP
}
// getPriorTheta
double getPriorTheta(double cv);
RcppExport SEXP _ggum_getPriorTheta(SEXP cvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type cv(cvSEXP);
    rcpp_result_gen = Rcpp::wrap(getPriorTheta(cv));
    return rcpp_result_gen;
END_RCPP
}
// getPriorAlpha
double getPriorAlpha(double cv);
RcppExport SEXP _ggum_getPriorAlpha(SEXP cvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type cv(cvSEXP);
    rcpp_result_gen = Rcpp::wrap(getPriorAlpha(cv));
    return rcpp_result_gen;
END_RCPP
}
// getPriorDelta
double getPriorDelta(double cv);
RcppExport SEXP _ggum_getPriorDelta(SEXP cvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type cv(cvSEXP);
    rcpp_result_gen = Rcpp::wrap(getPriorDelta(cv));
    return rcpp_result_gen;
END_RCPP
}
// getPriorTaus
double getPriorTaus(double cv);
RcppExport SEXP _ggum_getPriorTaus(SEXP cvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type cv(cvSEXP);
    rcpp_result_gen = Rcpp::wrap(getPriorTaus(cv));
    return rcpp_result_gen;
END_RCPP
}
// loglikelihoodRow
double loglikelihoodRow(IntegerVector responses, double theta, NumericVector alphas, NumericVector deltas, List taus);
RcppExport SEXP _ggum_loglikelihoodRow(SEXP responsesSEXP, SEXP thetaSEXP, SEXP alphasSEXP, SEXP deltasSEXP, SEXP tausSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type responses(responsesSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type deltas(deltasSEXP);
    Rcpp::traits::input_parameter< List >::type taus(tausSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikelihoodRow(responses, theta, alphas, deltas, taus));
    return rcpp_result_gen;
END_RCPP
}
// loglikelihoodCol
double loglikelihoodCol(IntegerVector responses, NumericVector thetas, double alpha, double delta, NumericVector taus);
RcppExport SEXP _ggum_loglikelihoodCol(SEXP responsesSEXP, SEXP thetasSEXP, SEXP alphaSEXP, SEXP deltaSEXP, SEXP tausSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type responses(responsesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type taus(tausSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikelihoodCol(responses, thetas, alpha, delta, taus));
    return rcpp_result_gen;
END_RCPP
}
// ggumMC3
NumericMatrix ggumMC3(IntegerMatrix data, int iters, int r_one, int r_two, int N, int W, Nullable<NumericVector> Temps);
RcppExport SEXP _ggum_ggumMC3(SEXP dataSEXP, SEXP itersSEXP, SEXP r_oneSEXP, SEXP r_twoSEXP, SEXP NSEXP, SEXP WSEXP, SEXP TempsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type iters(itersSEXP);
    Rcpp::traits::input_parameter< int >::type r_one(r_oneSEXP);
    Rcpp::traits::input_parameter< int >::type r_two(r_twoSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type W(WSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type Temps(TempsSEXP);
    rcpp_result_gen = Rcpp::wrap(ggumMC3(data, iters, r_one, r_two, N, W, Temps));
    return rcpp_result_gen;
END_RCPP
}
// ggumMCMC
NumericMatrix ggumMCMC(IntegerMatrix responseMatrix, IntegerVector Kvector, int iterations, int low, int high);
RcppExport SEXP _ggum_ggumMCMC(SEXP responseMatrixSEXP, SEXP KvectorSEXP, SEXP iterationsSEXP, SEXP lowSEXP, SEXP highSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type responseMatrix(responseMatrixSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Kvector(KvectorSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type low(lowSEXP);
    Rcpp::traits::input_parameter< int >::type high(highSEXP);
    rcpp_result_gen = Rcpp::wrap(ggumMCMC(responseMatrix, Kvector, iterations, low, high));
    return rcpp_result_gen;
END_RCPP
}
// ggumProbability
NumericVector ggumProbability(NumericVector k, double theta, double alpha, double delta, NumericVector tau);
RcppExport SEXP _ggum_ggumProbability(SEXP kSEXP, SEXP thetaSEXP, SEXP alphaSEXP, SEXP deltaSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(ggumProbability(k, theta, alpha, delta, tau));
    return rcpp_result_gen;
END_RCPP
}
// prob
double prob(const int choice, const double th, const double a, const double d, const NumericVector& t);
RcppExport SEXP _ggum_prob(SEXP choiceSEXP, SEXP thSEXP, SEXP aSEXP, SEXP dSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type choice(choiceSEXP);
    Rcpp::traits::input_parameter< const double >::type th(thSEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type d(dSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(prob(choice, th, a, d, t));
    return rcpp_result_gen;
END_RCPP
}
// probCol
NumericVector probCol(const IntegerVector& choices, const NumericVector& thetas, const double a, const double d, const NumericVector& t);
RcppExport SEXP _ggum_probCol(SEXP choicesSEXP, SEXP thetasSEXP, SEXP aSEXP, SEXP dSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type choices(choicesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type d(dSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(probCol(choices, thetas, a, d, t));
    return rcpp_result_gen;
END_RCPP
}
// probRow
NumericVector probRow(const IntegerVector& choices, const double th, const NumericVector& a, const NumericVector& d, const List& t);
RcppExport SEXP _ggum_probRow(SEXP choicesSEXP, SEXP thSEXP, SEXP aSEXP, SEXP dSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type choices(choicesSEXP);
    Rcpp::traits::input_parameter< const double >::type th(thSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const List& >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(probRow(choices, th, a, d, t));
    return rcpp_result_gen;
END_RCPP
}
// dlst
NumericVector dlst(NumericVector x, double df, double mu, double sigma);
RcppExport SEXP _ggum_dlst(SEXP xSEXP, SEXP dfSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(dlst(x, df, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// plst
NumericVector plst(NumericVector q, double df, double mu, double sigma);
RcppExport SEXP _ggum_plst(SEXP qSEXP, SEXP dfSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(plst(q, df, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// qlst
NumericVector qlst(NumericVector p, double df, double mu, double sigma);
RcppExport SEXP _ggum_qlst(SEXP pSEXP, SEXP dfSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(qlst(p, df, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// rlst
NumericVector rlst(int n, double df, double mu, double sigma);
RcppExport SEXP _ggum_rlst(SEXP nSEXP, SEXP dfSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(rlst(n, df, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// d_lst
double d_lst(double x, double df, double mu, double sigma);
RcppExport SEXP _ggum_d_lst(SEXP xSEXP, SEXP dfSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(d_lst(x, df, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// p_lst
double p_lst(double q, double df, double mu, double sigma);
RcppExport SEXP _ggum_p_lst(SEXP qSEXP, SEXP dfSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(p_lst(q, df, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// r_lst
double r_lst(double df, double mu, double sigma);
RcppExport SEXP _ggum_r_lst(SEXP dfSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(r_lst(df, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// update_theta_neg_MCMC
double update_theta_neg_MCMC(IntegerVector responses, double cv, NumericVector alphas, NumericVector deltas, List taus, double SD);
RcppExport SEXP _ggum_update_theta_neg_MCMC(SEXP responsesSEXP, SEXP cvSEXP, SEXP alphasSEXP, SEXP deltasSEXP, SEXP tausSEXP, SEXP SDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type responses(responsesSEXP);
    Rcpp::traits::input_parameter< double >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type deltas(deltasSEXP);
    Rcpp::traits::input_parameter< List >::type taus(tausSEXP);
    Rcpp::traits::input_parameter< double >::type SD(SDSEXP);
    rcpp_result_gen = Rcpp::wrap(update_theta_neg_MCMC(responses, cv, alphas, deltas, taus, SD));
    return rcpp_result_gen;
END_RCPP
}
// update_theta_pos_MCMC
double update_theta_pos_MCMC(IntegerVector responses, double cv, NumericVector alphas, NumericVector deltas, List taus, double SD);
RcppExport SEXP _ggum_update_theta_pos_MCMC(SEXP responsesSEXP, SEXP cvSEXP, SEXP alphasSEXP, SEXP deltasSEXP, SEXP tausSEXP, SEXP SDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type responses(responsesSEXP);
    Rcpp::traits::input_parameter< double >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type deltas(deltasSEXP);
    Rcpp::traits::input_parameter< List >::type taus(tausSEXP);
    Rcpp::traits::input_parameter< double >::type SD(SDSEXP);
    rcpp_result_gen = Rcpp::wrap(update_theta_pos_MCMC(responses, cv, alphas, deltas, taus, SD));
    return rcpp_result_gen;
END_RCPP
}
// update_alpha_MCMC
double update_alpha_MCMC(const IntegerVector& responses, const NumericVector& thetas, const double cv, const double delta, const NumericVector& taus, const double SD);
RcppExport SEXP _ggum_update_alpha_MCMC(SEXP responsesSEXP, SEXP thetasSEXP, SEXP cvSEXP, SEXP deltaSEXP, SEXP tausSEXP, SEXP SDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type responses(responsesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< const double >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< const double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type taus(tausSEXP);
    Rcpp::traits::input_parameter< const double >::type SD(SDSEXP);
    rcpp_result_gen = Rcpp::wrap(update_alpha_MCMC(responses, thetas, cv, delta, taus, SD));
    return rcpp_result_gen;
END_RCPP
}
// update_delta_MCMC
double update_delta_MCMC(const IntegerVector& responses, const NumericVector& thetas, const double alpha, const double cv, const NumericVector& taus, const double SD);
RcppExport SEXP _ggum_update_delta_MCMC(SEXP responsesSEXP, SEXP thetasSEXP, SEXP alphaSEXP, SEXP cvSEXP, SEXP tausSEXP, SEXP SDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type responses(responsesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type taus(tausSEXP);
    Rcpp::traits::input_parameter< const double >::type SD(SDSEXP);
    rcpp_result_gen = Rcpp::wrap(update_delta_MCMC(responses, thetas, alpha, cv, taus, SD));
    return rcpp_result_gen;
END_RCPP
}
// update_tau_MCMC
double update_tau_MCMC(const int k, const IntegerVector& responses, const NumericVector& thetas, const double alpha, const double delta, const NumericVector& taus, const double SD);
RcppExport SEXP _ggum_update_tau_MCMC(SEXP kSEXP, SEXP responsesSEXP, SEXP thetasSEXP, SEXP alphaSEXP, SEXP deltaSEXP, SEXP tausSEXP, SEXP SDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type responses(responsesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type taus(tausSEXP);
    Rcpp::traits::input_parameter< const double >::type SD(SDSEXP);
    rcpp_result_gen = Rcpp::wrap(update_tau_MCMC(k, responses, thetas, alpha, delta, taus, SD));
    return rcpp_result_gen;
END_RCPP
}
// update_theta_MC3
double update_theta_MC3(const double cv, const IntegerVector& choices, const NumericVector& a, const NumericVector& d, const List& t, const double temp, const double SD);
RcppExport SEXP _ggum_update_theta_MC3(SEXP cvSEXP, SEXP choicesSEXP, SEXP aSEXP, SEXP dSEXP, SEXP tSEXP, SEXP tempSEXP, SEXP SDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type choices(choicesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const List& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< const double >::type SD(SDSEXP);
    rcpp_result_gen = Rcpp::wrap(update_theta_MC3(cv, choices, a, d, t, temp, SD));
    return rcpp_result_gen;
END_RCPP
}
// update_alpha_MC3
double update_alpha_MC3(const double cv, const IntegerVector& choices, const NumericVector& th, const double d, const NumericVector& t, const double temp, const double SD);
RcppExport SEXP _ggum_update_alpha_MC3(SEXP cvSEXP, SEXP choicesSEXP, SEXP thSEXP, SEXP dSEXP, SEXP tSEXP, SEXP tempSEXP, SEXP SDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type choices(choicesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type th(thSEXP);
    Rcpp::traits::input_parameter< const double >::type d(dSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< const double >::type SD(SDSEXP);
    rcpp_result_gen = Rcpp::wrap(update_alpha_MC3(cv, choices, th, d, t, temp, SD));
    return rcpp_result_gen;
END_RCPP
}
// update_delta_MC3
double update_delta_MC3(const double cv, const IntegerVector& choices, const NumericVector& th, const double a, const NumericVector& t, const double temp, const double SD);
RcppExport SEXP _ggum_update_delta_MC3(SEXP cvSEXP, SEXP choicesSEXP, SEXP thSEXP, SEXP aSEXP, SEXP tSEXP, SEXP tempSEXP, SEXP SDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type choices(choicesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type th(thSEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< const double >::type SD(SDSEXP);
    rcpp_result_gen = Rcpp::wrap(update_delta_MC3(cv, choices, th, a, t, temp, SD));
    return rcpp_result_gen;
END_RCPP
}
// update_tau_MC3
double update_tau_MC3(const int k, const IntegerVector& choices, const NumericVector& th, const double a, const double d, const NumericVector& t, const double temp, const double SD);
RcppExport SEXP _ggum_update_tau_MC3(SEXP kSEXP, SEXP choicesSEXP, SEXP thSEXP, SEXP aSEXP, SEXP dSEXP, SEXP tSEXP, SEXP tempSEXP, SEXP SDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type choices(choicesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type th(thSEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type d(dSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< const double >::type SD(SDSEXP);
    rcpp_result_gen = Rcpp::wrap(update_tau_MC3(k, choices, th, a, d, t, temp, SD));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ggum_d4beta", (DL_FUNC) &_ggum_d4beta, 5},
    {"_ggum_p4beta", (DL_FUNC) &_ggum_p4beta, 5},
    {"_ggum_q4beta", (DL_FUNC) &_ggum_q4beta, 5},
    {"_ggum_r4beta", (DL_FUNC) &_ggum_r4beta, 5},
    {"_ggum_d_4beta", (DL_FUNC) &_ggum_d_4beta, 5},
    {"_ggum_r_4beta", (DL_FUNC) &_ggum_r_4beta, 4},
    {"_ggum_getPriorTheta", (DL_FUNC) &_ggum_getPriorTheta, 1},
    {"_ggum_getPriorAlpha", (DL_FUNC) &_ggum_getPriorAlpha, 1},
    {"_ggum_getPriorDelta", (DL_FUNC) &_ggum_getPriorDelta, 1},
    {"_ggum_getPriorTaus", (DL_FUNC) &_ggum_getPriorTaus, 1},
    {"_ggum_loglikelihoodRow", (DL_FUNC) &_ggum_loglikelihoodRow, 5},
    {"_ggum_loglikelihoodCol", (DL_FUNC) &_ggum_loglikelihoodCol, 5},
    {"_ggum_ggumMC3", (DL_FUNC) &_ggum_ggumMC3, 7},
    {"_ggum_ggumMCMC", (DL_FUNC) &_ggum_ggumMCMC, 5},
    {"_ggum_ggumProbability", (DL_FUNC) &_ggum_ggumProbability, 5},
    {"_ggum_prob", (DL_FUNC) &_ggum_prob, 5},
    {"_ggum_probCol", (DL_FUNC) &_ggum_probCol, 5},
    {"_ggum_probRow", (DL_FUNC) &_ggum_probRow, 5},
    {"_ggum_dlst", (DL_FUNC) &_ggum_dlst, 4},
    {"_ggum_plst", (DL_FUNC) &_ggum_plst, 4},
    {"_ggum_qlst", (DL_FUNC) &_ggum_qlst, 4},
    {"_ggum_rlst", (DL_FUNC) &_ggum_rlst, 4},
    {"_ggum_d_lst", (DL_FUNC) &_ggum_d_lst, 4},
    {"_ggum_p_lst", (DL_FUNC) &_ggum_p_lst, 4},
    {"_ggum_r_lst", (DL_FUNC) &_ggum_r_lst, 3},
    {"_ggum_update_theta_neg_MCMC", (DL_FUNC) &_ggum_update_theta_neg_MCMC, 6},
    {"_ggum_update_theta_pos_MCMC", (DL_FUNC) &_ggum_update_theta_pos_MCMC, 6},
    {"_ggum_update_alpha_MCMC", (DL_FUNC) &_ggum_update_alpha_MCMC, 6},
    {"_ggum_update_delta_MCMC", (DL_FUNC) &_ggum_update_delta_MCMC, 6},
    {"_ggum_update_tau_MCMC", (DL_FUNC) &_ggum_update_tau_MCMC, 7},
    {"_ggum_update_theta_MC3", (DL_FUNC) &_ggum_update_theta_MC3, 7},
    {"_ggum_update_alpha_MC3", (DL_FUNC) &_ggum_update_alpha_MC3, 7},
    {"_ggum_update_delta_MC3", (DL_FUNC) &_ggum_update_delta_MC3, 7},
    {"_ggum_update_tau_MC3", (DL_FUNC) &_ggum_update_tau_MC3, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_ggum(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
