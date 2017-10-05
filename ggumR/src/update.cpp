#include "../inst/include/ggumR.h"

using namespace Rcpp;

//[[Rcpp::export]]
double dhnormpos(double x){
    if ( x < 0 ) {
        return 0.0;
    }
    return 2.0 * R::dnorm(x, 0.0, 1.0, 0);
}


//[[Rcpp::export]]
double dhnormneg(double x){
    if ( x > 0 ) {
        return 0.0;
    }
    return 2.0 * R::dnorm(x, 0.0, 1.0, 0);
}


//[[Rcpp::export]]
double updateTheta(double cv, IntegerVector choices, NumericVector a,
        NumericVector d, List t, double temp){
    double pv = r_lst(1, cv, 1);
    double cvPrior = R::dnorm(cv, 0.0, 1.0, 0);
    double pvPrior = R::dnorm(pv, 0.0, 1.0, 0);
    double cvL = sum(log(probRow(choices, cv, a, d, t)));
    double pvL = sum(log(probRow(choices, pv, a, d, t)));
    double r = pow((exp(pvL - cvL)) * (pvPrior/cvPrior), temp);
    if ( r > 1 || R::runif(0, 1) < r) {
        return pv;
    }
    return cv;
}


//[[Rcpp::export]]
double updateThetaPos(double cv, IntegerVector choices, NumericVector a,
        NumericVector d, List t, double temp){
    double pv = r_lst(1, cv, 1);
    double cvPrior = dhnormpos(cv);
    double pvPrior = dhnormpos(pv);
    double cvL = sum(log(probRow(choices, cv, a, d, t)));
    double pvL = sum(log(probRow(choices, pv, a, d, t)));
    double r = pow((exp(pvL - cvL)) * (pvPrior/cvPrior), temp);
    if ( r > 1 || R::runif(0, 1) < r) {
        return pv;
    }
    return cv;
}


//[[Rcpp::export]]
double updateThetaNeg(double cv, IntegerVector choices, NumericVector a,
        NumericVector d, List t, double temp){
    double pv = r_lst(1, cv, 1);
    double cvPrior = dhnormneg(cv);
    double pvPrior = dhnormneg(pv);
    double cvL = sum(log(probRow(choices, cv, a, d, t)));
    double pvL = sum(log(probRow(choices, pv, a, d, t)));
    double r = pow((exp(pvL - cvL)) * (pvPrior/cvPrior), temp);
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
    double r = pow((exp(pvL - cvL)) * (pvPrior/cvPrior), temp);
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
    double r = pow((exp(pvL - cvL)) * (pvPrior/cvPrior), temp);
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
    double r = pow((exp(pvL - cvL)) * (pvPrior/cvPrior), temp);
    if ( r > 1 || R::runif(0, 1) < r) {
        return pv[k];
    }
    return t[k];
}

