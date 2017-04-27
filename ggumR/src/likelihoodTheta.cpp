#include "../inst/include/ggumR.h"

using namespace Rcpp;

//' GGUM LogLikelihood - Theta
//' 
//' Calculate the Log Likelihood parameters Theta for the GGUM.
//' 
//' @param theta A numeric vector of length one giving the individual's latent
//'   trait parameter
//' @param responseVector A numeric vector of length N (the number of
//'   respondents) giving the option chosen the individual for each item j
//' @param alphas A numeric vector of length n; each element of the vector is an
//'   item's discrimination parameter
//' @param deltas A numeric vector of length n; each element of the vector is an
//'   item's location parameter
//' @param taus A list of numeric vectors; each list element j is a numeric
//'   vector of threshold parameters for item j's options (where the first
//'   element of the vector should be zero).
//' @return The (log) likelihood of the vector of interest.
//' @rdname ggumLogLikelihoodTheta
//' @export
//[[Rcpp::export]]
double loglikelihoodTheta(NumericVector responses, double theta,
                       NumericVector alphas, NumericVector deltas, List taus){
   NumericVector thetaProbs = probRow(responses, theta, alphas, deltas, taus);
   double llTheta = sum(log(thetaProbs));
   return llTheta;
}
