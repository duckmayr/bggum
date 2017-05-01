#include "../inst/include/ggumR.h"

using namespace Rcpp;

//' GGUM Log Likelihood
//' 
//' Calculate the log likelihood of data for the GGUM given parameter values.
//'
//' This function calculates the log likelihood of a \bold{vector} of
//' responses given values for the parameters relevant to the responses.
//' This could be all of a respondent \eqn{i}'s responses to every item \eqn{j},
//' or all respondents' responses to an item \eqn{j}. We calculate likelihood
//' of vectors rather than the likelihood of the entire response matrix since
//' when computing acceptance ratios for the MCMC algorithm, we divide products
//' of probabilities, most of which will cancel out -- thus, we only ever need
//' the likelihood of one vector at any given time.
//' 
//' @param theta For \code{loglikelihoodRow}, a numeric vector of length one
//'   giving the individual's latent trait parameter
//' @param responses For \code{loglikelihoodRow}, a numeric vector of length
//'   length n (the number of items) giving the option chosen by the individual
//'   for each item j; for \code{loglikelihoodCol}, a numeric vector of length
//'   N (the number of respondents), giving the option chosen by each individual
//'   i to the item
//' @param alphas A numeric vector of length n; each element of the vector is an
//'   item's discrimination parameter
//' @param deltas A numeric vector of length n; each element of the vector is an
//'   item's location parameter
//' @param taus For \code{loglikelihoodRow}, a list of numeric vectors where
//'   each list element j is a numeric vector of threshold parameters for item
//'   j's options (where the first element of the vector should be zero);
//'   for \code{loglikelihoodCol}, this is only the numeric vector of threshold
//'   parameters for the item of interest
//' @param thetas A numeric vector of length N, each each element of which is
//'   an individual's latent trait parameter
//' @param alpha A numeric vector of length one giving the item's
//'   discrimination parameter
//' @param delta A numeric vector of length one giving the item's
//'   location parameter
//'
//' @return The (log) likelihood of the vector of interest.
//' @rdname ggumLogLikelihood
//' @export
//[[Rcpp::export]]
double loglikelihoodRow(NumericVector responses, double theta,
        NumericVector alphas, NumericVector deltas, List taus){
    return sum(log(na_omit(probRow(responses, theta, alphas, deltas, taus))));
}

//' @rdname ggumLogLikelihood
//' @export
//[[Rcpp::export]]
double loglikelihoodCol(NumericVector responses, NumericVector thetas,
        double alpha, double delta, NumericVector taus){
    return sum(log(na_omit(probCol(responses, thetas, alpha, delta, taus))));
}

