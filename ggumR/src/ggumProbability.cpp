#include <Rcpp.h>

using namespace Rcpp;

double prob(int k, double theta, double alpha, double delta, 
            NumericVector tau){
    int K = tau.size();
    NumericVector cumSum = cumsum(tau);
    NumericVector numerator(K);
    for ( int x = 0; x < K; x++ ) {
        numerator[x] = exp(alpha * (x * (theta - delta) - cumSum[x]))
            + exp(alpha * ((2*K - 1 - x) * (theta - delta) - cumSum[x]));
    }
    double denominator = sum(numerator);
    return numerator[k-1]/denominator;
}

NumericVector probRow(NumericVector responseVec, double theta,
                      NumericVector alphas, NumericVector deltas, List taus){
    int m = responseVec.size();
    NumericVector result(m);
    for ( int j = 0; j < m; j++ ) {
        if ( ISNA(responseVec[j]) ) {
            result[j] = NA_REAL;
            continue;
        }
        NumericVector tau = as<NumericVector>(taus[j]);
        result[j] = prob(responseVec[j], theta, alphas[j], deltas[j], tau);
    }
    return result;
}

NumericVector probCol(NumericVector responseVec, NumericVector thetas,
                      double alpha, double delta, NumericVector taus){
    int n = responseVec.size();
    NumericVector result(n);
    for ( int i = 0; i < n; i++ ) {
        if ( ISNA(responseVec[i]) ) {
            result[i] = NA_REAL;
            continue;
        }
        result[i] = prob(responseVec[i], thetas[i], alpha, delta, taus);
    }
    return result;
}

//' GGUM Probability Function
//'
//' Calculate the probability of a response according to the GGUM
//' 
//' The General Graded Unfolding Model (GGUM) is an item response model
//' designed to consider the possibility of disagreement for opposite reasons.
//' This function gives the probability of a respondent's response to a test
//' item given item and respondent parameters. The user can calculate the
//' probability of one particular response or for any number of the possible
//' responses to the item (see the \code{k} parameter described below).
//' 
//' @param k A numeric vector giving the item(s) for which response
//'   probability should be calculated
//' @param theta A numeric vector of length one giving the respondent's latent
//'   trait score
//' @param alpha A numeric vector of length one giving the item's
//'   discrimination parameter
//' @param delta A numeric vector of length one giving the item's location
//'   parameter
//' @param tau A numeric vector of length K (the number of possible responses)
//'   giving the options' threshold parameters; the first element of \code{tau}
//'   should be zero
//'
//' @return A numeric vector the same length of \code{k} giving the response
//'   probabilities
//' @rdname ggumProbability
//' @export
//[[Rcpp::export]]
NumericVector ggumProbability(NumericVector k, double theta, double alpha,
        double delta, NumericVector tau){
    int K = tau.size();
    NumericVector cumSum = cumsum(tau);
    NumericVector result(K);
    for ( int x = 0; x < K; x++ ) {
        result[x] = exp(alpha * (x * (theta - delta) - cumSum[x]))
            + exp(alpha * ((2*K - 1 - x) * (theta - delta) - cumSum[x]));
    }
    result = result / sum(result);
    return result[k-1];
}

