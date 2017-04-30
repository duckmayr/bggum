#include "../inst/include/ggumR.h"

using namespace Rcpp;

//' GGUM acceptance
//'
//' Computes an acceptance ratio for the parameters used in the MCMC algorithm
//' for the GGUM.
//' 
//' @param responses A numeric vector giving the responses by the respondent
//'   or item of interest; for \code{acceptanceTheta} this would be one row of
//'   the response matrix, while for \code{acceptanceAlpha},
//'   \code{acceptanceDelta}, or \code{acceptanceTau} this would be one column
//'   of the response matrix
//' @param cv The current value of the paramenter of interest
//' @param thetas The numeric vector of current values of the items' theta
//'   parameters
//' @param alphas The numeric vector of current values of the items' alpha
//'   parameters
//' @param alpha The alpha parameter for the item of interest (for
//'   \code{acceptanceDelta} and \code{acceptanceTau}, there is only one
//'   alpha needed for the calculations)
//' @param deltas The numeric vector of current values of the items' delta
//'   parameters
//' @param delta The delta parameter for the item of interest (for
//'   \code{acceptanceAlpha} and \code{acceptanceTau}, there is only one
//'   delta needed for the calculations)
//' @param taus The list or vector of current values of the items' tau
//'   parameters (for \code{acceptanceTheta} this is a list because there are
//'   multiple tau vectors that must be used, but for the other acceptance
//'   functions, only one numeric vector is needed)
//' @param SD The sigma parameter to be used for the proposal distribution
//'
//' @return If the acceptance ratio is greater than one or greater than a
//'   number randomly generated from the standard uniform distribution,
//'   the proposed value is returned; otherwise, the \code{cv} is returned.
//' @export
//[[Rcpp::export]]
double acceptanceTheta(NumericVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD){
    // pv stands for proposed value, cv for current value
    // First we propose a new value for the current theta;
    // we use a truncated location-scale T distribution:
    double pv = r_trunclst(1, cv, SD, -10, 10);
    // Then we get the prior probabilities that theta is equal to
    // the current value (cvPrior) and proposed value (pvPrior):
    double pvPrior = R::dnorm(pv, 0, 1, 0);
    double cvPrior = R::dnorm(cv, 0, 1, 0);
    // And the vector of probabilities that we would observe the values in
    // responses given the current value of theta and the propposed value:
    NumericVector cvProbs = probRow(responses, cv, alphas, deltas, taus);
    NumericVector pvProbs = probRow(responses, pv, alphas, deltas, taus);
    // Then we calculate the acceptance rate:
    double acceptRate = exp(sum(log(pvProbs)) + log(pvPrior)
            - sum(log(cvProbs)) - log(cvPrior));
    // If the acceptance rate is greater than one, not a number,
    // or larger than a uniform deviate, we accept the proposal:
    if ( acceptRate > 1 || isnan(acceptRate) || R::runif(0, 1) < acceptRate) {
        return pv;
    }
    // Otherwise we reject the proposal
    return cv;
}

//' @export
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

//' @export
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

//' @export
//[[Rcpp::export]]
double acceptanceTau(int k, NumericVector responses, NumericVector thetas,
        double alpha, double delta, NumericVector taus, double SD){
    // For taus, we need a copy of the entire vector
    NumericVector pv = clone(taus);
    // So that when we replace one of its values with a proposed value,
    // we can still compute probability of responses according to the GGUM
    // under the original tau vector
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

