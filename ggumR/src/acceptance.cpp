#include "../inst/include/ggumR.h"

using namespace Rcpp;

//' GGUM acceptance
//'
//' Determines a new value for the parameters used in the MCMC algorithm
//' for the GGUM.
//' 
//' Given a current value \code{cv} for the parameter of interest,
//' as well as the responses and other parameters relevant to estimating the
//' parameter of interest, a new proposal for the parameter is generated
//' and accepted with probability
//' \deqn{\min \{1, \frac{\mathcal{L}(X|\theta^*)\pi(\theta^*)}{%
//' \mathcal{L}(X|\theta)\pi(\theta)}\}}{%
//' min {1, (L(X|\theta*)\pi(\theta*) / (L(X|\theta)\pi(\theta))}}
//' where \eqn{\theta^*}{\theta*} is the proposed value,
//' \eqn{\theta} is the current value, \eqn{\pi(\cdot)}{\pi(.)}
//' is the prior probability a parameter takes a value, and
//' \eqn{\mathcal{L}(X|\cdot)}{L(X|.)} is the likelihood of observing the
//' responses given a parameter value.
//'
//' @section Warning:
//' For \code{acceptanceTau}, the tau vector is indexed beginning with 0!
//' When using this function, you must supply an index of one less than the
//' index you would supply if subsetting in \code{R}.
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
//' @param k For \code{acceptanceTau}, an integer giving the index
//'   (NOTE: the index must be specified as the vector would be indexed in
//'   C++; i.e., the R index - 1 since C++ begins indexing at 0) in the
//'   tau vector of the tau parameter of interest.
//'
//' @return If the acceptance ratio is greater than one or greater than a
//'   number randomly generated from the standard uniform distribution,
//'   the proposed value is returned; otherwise, the \code{cv} is returned.
//'
//' @name GGUM MCMC Proposal Acceptance
//' @aliases acceptanceTheta, acceptanceAlpha, acceptanceDelta, acceptanceTau,
//'   acceptance
//' @rdname acceptance
//' @export
//[[Rcpp::export]]
double acceptanceTheta(IntegerVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD){
    // pv stands for proposed value, cv for current value
    // First we propose a new value for the current theta;
    // we use a truncated location-scale T distribution:
    double pv = r_lst(1, cv, SD);
    // Then we get the prior probabilities that theta is equal to
    // the current value (cvPrior) and proposed value (pvPrior):
    double pvPrior = R::dnorm(pv, 0, 1, 0);
    double cvPrior = R::dnorm(cv, 0, 1, 0);
    // The (log) likelihood of the data given the current and proposed thetas:
    double cvL = sum(log(probRow(responses, cv, alphas, deltas, taus)));
    double pvL = sum(log(probRow(responses, pv, alphas, deltas, taus)));
    // Then we calculate the acceptance rate:
    double acceptRate = exp(pvL - cvL) * (pvPrior/cvPrior);
    // If the acceptance rate is greater than one, not a number,
    // or larger than a uniform deviate, we accept the proposal:
    if ( acceptRate > 1 || R::runif(0, 1) < acceptRate) {
        return pv;
    }
    // Otherwise we reject the proposal
    return cv;
}

//[[Rcpp::export]]
double acceptanceThetaNeg(IntegerVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD){
    double pv = r_lst(1, cv, SD);
    double pvPrior = R::dnorm(pv, 0, 1, 0);
    double cvPrior = R::dnorm(cv, 0, 1, 0);
    double cvL = sum(log(probRow(responses, cv, alphas, deltas, taus)));
    double pvL = sum(log(probRow(responses, pv, alphas, deltas, taus)));
    double acceptRate = exp(pvL - cvL) * (pvPrior/cvPrior);
    if ( acceptRate > 1 || R::runif(0, 1) < acceptRate) {
        return pv;
    }
    return cv;
}

//[[Rcpp::export]]
double acceptanceThetaPos(IntegerVector responses, double cv,
        NumericVector alphas, NumericVector deltas, List taus, double SD){
    double pv = r_lst(1, cv, SD);
    double pvPrior = R::dnorm(pv, 0, 1, 0);
    double cvPrior = R::dnorm(cv, 0, 1, 0);
    double cvL = sum(log(probRow(responses, cv, alphas, deltas, taus)));
    double pvL = sum(log(probRow(responses, pv, alphas, deltas, taus)));
    double acceptRate = exp(pvL - cvL) * (pvPrior/cvPrior);
    if ( acceptRate > 1 || R::runif(0, 1) < acceptRate) {
        return pv;
    }
    return cv;
}

//' @rdname acceptance
//' @export
//[[Rcpp::export]]
double acceptanceAlpha(IntegerVector responses, NumericVector thetas,
        double cv, double delta, NumericVector taus, double SD){
    double pv = r_lst(1, cv, SD);
    double pvPrior = d_4beta(pv, 1.5, 1.5, 0.25, 4);
    double cvPrior = d_4beta(cv, 1.5, 1.5, 0.25, 4);
    double cvL = sum(log(probCol(responses, thetas, cv, delta, taus)));
    double pvL = sum(log(probCol(responses, thetas, pv, delta, taus)));
    double acceptRate = exp(pvL - cvL) * (pvPrior/cvPrior);
    if ( acceptRate > 1 || R::runif(0, 1) < acceptRate) {
        return pv;
    }
    return cv;
}

//' @rdname acceptance
//' @export
//[[Rcpp::export]]
double acceptanceDelta(IntegerVector responses, NumericVector thetas,
        double alpha, double cv, NumericVector taus, double SD){
    double pv = r_lst(1, cv, SD);
    double pvPrior = d_4beta(pv, 2, 2, -5, 5);
    double cvPrior = d_4beta(cv, 2, 2, -5, 5);
    double cvL = sum(log(probCol(responses, thetas, alpha, cv, taus)));
    double pvL = sum(log(probCol(responses, thetas, alpha, pv, taus)));
    double acceptRate = exp(pvL - cvL) * (pvPrior/cvPrior);
    if ( acceptRate > 1 || R::runif(0, 1) < acceptRate) {
        return pv;
    }
    return cv;
}

//' @rdname acceptance
//' @export
//[[Rcpp::export]]
double acceptanceTau(int k, IntegerVector responses, NumericVector thetas,
        double alpha, double delta, NumericVector taus, double SD){
    // For taus, we need a copy of the entire vector
    NumericVector pv = clone(taus);
    // So that when we replace one of its values with a proposed value,
    // we can still compute probability of responses according to the GGUM
    // under the original tau vector
    pv[k] = r_lst(1, taus[k], SD);
    double pvPrior = d_4beta(pv[k], 2, 2, -2, 0);
    double cvPrior = d_4beta(taus[k], 2, 2, -2, 0);
    double cvL = sum(log(probCol(responses, thetas, alpha, delta, taus)));
    double pvL = sum(log(probCol(responses, thetas, alpha, delta, pv)));
    double acceptRate = exp(pvL - cvL) * (pvPrior/cvPrior);
    if ( acceptRate > 1 || R::runif(0, 1) < acceptRate) {
        return pv[k];
    }
    return taus[k];
}

