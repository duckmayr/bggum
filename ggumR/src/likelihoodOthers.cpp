#include "../inst/include/ggumR.h"

using namespace Rcpp;

//' GGUM LogLikelihood Alpha, Delta, and Taus.
//' 
//' Calculate the Log Likelihood parameters Alpha, Delta, and Taus
//' for the GGUM.
//' 
//' @param thetas A numeric vector of length N (the number of respondents); each
//'   each element of the vector is an individual's latent trait parameter
//' @param responseVector A numeric vector of length n (the number of items)
//'   giving the option chosen by individual i for the item
//' @param alpha A numeric vector of length one giving the item's discrimination
//'   parameter
//' @param delta A numeric vector of length one giving the item's location
//'   parameter
//' @param taus A numeric vector of length K (the number of options) giving the
//'   threshold parameter for each option; the first element should be zero.
//' 
//' @return The (log) likelihood of the vector of interest.
//' @rdname ggumLogLikelihoodOthers
//' @export
//[[Rcpp::export]]
double loglikelihoodAlpha(NumericVector responses, NumericVector thetas,
                          double alpha, double delta, NumericVector taus){
   NumericVector AlphaProbs = probCol(responses, thetas, alpha, delta, taus);
   double llAlpha = sum(log(AlphaProbs));
   return llAlpha;
}

//[[Rcpp::export]]
double loglikelihoodDelta(NumericVector responses, NumericVector thetas,
                          double alpha, double delta, NumericVector taus){
   NumericVector DeltaProbs = probCol(responses, thetas, alpha, delta, taus);
   double llDelta = sum(log(DeltaProbs));
   return llDelta;
}

//[[Rcpp::export]]
double loglikelihoodTau(NumericVector responses, NumericVector thetas,
                        double alpha, double delta, NumericVector taus){
   NumericVector TausProbs = probCol(responses, thetas, alpha, delta, taus);
   double llTaus = sum(log(TausProbs));
   return llTaus;
}