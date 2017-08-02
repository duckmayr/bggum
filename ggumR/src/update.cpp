#include "../inst/include/ggumR.h"

using namespace Rcpp;

//[[Rcpp::export]]
double updateTheta(double cv, IntegerVector choices, NumericVector a,
        NumericVector d, List t, double temp, double hi, double lo){
    double pv = r_lst(1, cv, 1);
    double cvPrior, pvPrior;
    double scale = R::pnorm(hi, 0, 1, 1, 0) - R::pnorm(lo, 0, 1, 1, 0);
    if ( cv < lo || cv > hi ) {
        return pv;
    }
    else {
        cvPrior = R::dnorm(cv, 0, 1, 0) / scale;
    }
    if ( pv < lo || pv > hi ) {
        pvPrior = 0;
    }
    else {
        pvPrior = R::dnorm(pv, 0, 1, 0) / scale;
    }
    double cvL = sum(log(probRow(choices, cv, a, d, t)));
    double pvL = sum(log(probRow(choices, pv, a, d, t)));
    double r = pow((exp(pvL - cvL)) * (pvPrior/cvPrior), 1/temp);
    if ( r > 1 || R::runif(0, 1) < r) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double updateAlpha(double cv, IntegerVector choices, NumericVector th,
        double d, NumericVector t, double temp){
    double pv = r_lst(1, cv, 0.15);
    double pvPrior = d_4beta(pv, 1.5, 1.5, 0.25, 4);
    double cvPrior = d_4beta(cv, 1.5, 1.5, 0.25, 4);
    double cvL = sum(log(probCol(choices, th, cv, d, t)));
    double pvL = sum(log(probCol(choices, th, pv, d, t)));
    double r = pow((exp(pvL - cvL)) * (pvPrior/cvPrior), 1/temp);
    if ( r > 1 || R::runif(0, 1) < r) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double updateDelta(double cv, IntegerVector choices, NumericVector th,
        double a, NumericVector t, double temp){
    double pv = r_lst(1, cv, 0.2);
    double pvPrior = d_4beta(pv, 2, 2, -5, 5);
    double cvPrior = d_4beta(cv, 2, 2, -5, 5);
    double cvL = sum(log(probCol(choices, th, a, cv, t)));
    double pvL = sum(log(probCol(choices, th, a, pv, t)));
    double r = pow((exp(pvL - cvL)) * (pvPrior/cvPrior), 1/temp);
    if ( r > 1 || R::runif(0, 1) < r) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double updateTau(int k, IntegerVector choices, NumericVector th,
        double a, double d, NumericVector t, double temp){
    NumericVector pv = clone(t);
    pv[k] = r_lst(1, t[k], 0.5);
    double pvPrior = d_4beta(pv[k], 2, 2, -6, 6);
    double cvPrior = d_4beta(t[k], 2, 2, -6, 6);
    double cvL = sum(log(probCol(choices, th, a, d, t)));
    double pvL = sum(log(probCol(choices, th, a, d, pv)));
    double r = pow((exp(pvL - cvL)) * (pvPrior/cvPrior), 1/temp);
    if ( r > 1 || R::runif(0, 1) < r) {
        return pv[k];
    }
    return t[k];
}

