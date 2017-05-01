#include <Rcpp.h>

using namespace Rcpp;

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
//' The probability that respondent \eqn{i} chooses option \eqn{k} for item
//' \eqn{j} is given by
//' \deqn{\frac{\exp (\alpha_j [k (\theta_i - \delta_j) -
//' \sum_{m=0}^k \tau_{jm}]) + \exp (\alpha_j [(2K - k - 1)
//' (\theta_i - \delta_j) - \sum_{m=0}^k \tau_{jm}])}{%
//' \sum_{l=0}^{K-1} [\exp (\alpha_j [l (\theta_i - \delta_j) -
//' \sum_{m=0}^l \tau_{jm}]) + \exp (\alpha_j [(2K - l - 1)
//' (\theta_i - \delta_j) - \sum_{m=0}^l \tau_{jm}])]}}{%
//' (exp(\alpha_j [k(\theta_i-\delta_j) - 
//' \Sigma_{m=0}^k \tau_{jm}]) + exp(\alpha_j [(2K - k - 1)
//' (\theta_i - \delta_j) - \Sigma_{m=0}^k \tau_{jm}])) /
//' (\Sigma_{l=0}^{K-1} [exp (\alpha_j [l (\theta_i - \delta_j) -
//' \Sigma_{m=0}^l \tau_{jm}]) + exp (\alpha_j [(2K - l - 1)
//' (\theta_i - \delta_j) - \Sigma_{m=0}^l \tau_{jm}])])},
//' where \eqn{\theta_i} is \eqn{i}'s latent trait parameter,
//' \eqn{\alpha_j} is the item's discrimination paramter,
//' \eqn{\delta_j} is the item's location paramter,
//' \eqn{\tau_{j0}, \ldots, \tau_{j(K-1)}} are the options' threshold
//' parameters, and \eqn{\tau_{j0}} is 0,
//' \eqn{K} is the number of options for item \eqn{j}, and
//' the options are indexed by \eqn{k = 0, \ldots, K-1}.
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
//'
//' @references Roberts, James S., John R. Donoghue, and James E. Laughlin.
//'   2000. ``A General Item Response Theory Model for Unfolding
//'   Unidimensional Polytomous Responses." \emph{Applied Psychological
//'   Measurement} 24(1): 3--32.
//' @references de la Torre, Jimmy, Stephen Stark, and Oleksandr S.
//'   Chernyshenko. 2006. ``Markov Chain Monte Carlo Estimation of Item
//'   Parameters for the Generalized Graded Unfolding Model." \emph{Applied
//'   Psychological Measurement} 30(3): 216--232.
//'
//' @rdname ggumProbability
//' @export
//[[Rcpp::export]]
NumericVector ggumProbability(NumericVector k, double theta, double alpha,
        double delta, NumericVector tau){
    // K is the number of options
    int K = tau.size();
    // The elements of the cumulative summation of the tau vector is used
    // multiple times in calculating response probabilities,
    // so we compute and store it now
    NumericVector cumSum = cumsum(tau);
    // Then we create a numeric vector to calculate the probability of
    // every possible response; this is done because the denominator of
    // the probability calculation is the sum of every possible numerator
    // value, so even if the user does not want the probability of all
    // possible responses, we must compute these values anyway
    NumericVector result(K);
    // For each element of that vector, we compute the numerator expression
    for ( int x = 0; x < K; x++ ) {
        result[x] = exp(alpha * (x * (theta - delta) - cumSum[x]))
            + exp(alpha * ((2*K - 1 - x) * (theta - delta) - cumSum[x]));
    }
    // Then we divide by the sum
    result = result / sum(result);
    // and return the kth element(s) (k-1 to account for R-C++ differences)
    return result[k-1];
}


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
