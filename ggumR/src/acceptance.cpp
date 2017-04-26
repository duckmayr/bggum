#include "../inst/include/ggumR.h"

using namespace Rcpp;

//[[Rcpp::export]]
double acceptanceTheta(NumericVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD){
    // pv stands for proposed value, cv for current value
    double pv = r_trunclst(1, cv, SD, -10, 10);
    double pvPrior = R::dnorm(pv, 0, 1, 0);
    double cvPrior = R::dnorm(cv, 0, 1, 0);
    NumericVector cvProbs = probRow(responses, cv, alphas, deltas, taus);
    NumericVector pvProbs = probRow(responses, pv, alphas, deltas, taus);
    double acceptRate = exp(sum(log(pvProbs)) + log(pvPrior)
            - sum(log(cvProbs)) - log(cvPrior));
    if ( acceptRate > 1 || isnan(acceptRate) || R::runif(0, 1) < acceptRate) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double acceptanceAlpha(NumericVector responses, NumericVector thetas,
        double cv, double delta, NumericVector taus, double SD){
    double pv = r_trunclst(1, cv, SD, 0.25, 4);
    double pvPrior = d_4beta(pv, 1.5, 1.5, 0.25, 4);
    double cvPrior = d_4beta(cv, 1.5, 1.5, 0.25, 4);
    NumericVector cvProbs = probCol(responses, thetas, cv, delta, taus);
    NumericVector pvProbs = probCol(responses, thetas, pv, delta, taus);
    double acceptRate = exp(sum(log(pvProbs)) + log(pvPrior)
            - sum(log(cvProbs)) - log(cvPrior));
    if ( acceptRate > 1 || isnan(acceptRate) || R::runif(0, 1) < acceptRate) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double acceptanceDelta(NumericVector responses, NumericVector thetas,
        double alpha, double cv, NumericVector taus, double SD){
    double pv = r_trunclst(1, cv, SD, -5, 5);
    double pvPrior = d_4beta(pv, 2, 2, -5, 5);
    double cvPrior = d_4beta(cv, 2, 2, -5, 5);
    NumericVector cvProbs = probCol(responses, thetas, alpha, cv, taus);
    NumericVector pvProbs = probCol(responses, thetas, alpha, pv, taus);
    double acceptRate = exp(sum(log(pvProbs)) + log(pvPrior)
            - sum(log(cvProbs)) - log(cvPrior));
    if ( acceptRate > 1 || isnan(acceptRate) || R::runif(0, 1) < acceptRate) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double acceptanceTau(int k, NumericVector responses, NumericVector thetas,
        double alpha, double delta, NumericVector taus, double SD){
    NumericVector pv = clone(taus);
    pv[k] = r_trunclst(1, taus[k], SD, -6, 6);
    double pvPrior = d_4beta(pv[k], 2, 2, -6, 6);
    double cvPrior = d_4beta(taus[k], 2, 2, -6, 6);
    NumericVector cvProbs = probCol(responses, thetas, alpha, delta, taus);
    NumericVector pvProbs = probCol(responses, thetas, alpha, delta, pv);
    double acceptRate = exp(sum(log(pvProbs)) + log(pvPrior)
            - sum(log(cvProbs)) - log(cvPrior));
    if ( acceptRate > 1 || isnan(acceptRate) || R::runif(0, 1) < acceptRate) {
        return pv[k];
    }
    return taus[k];
}

