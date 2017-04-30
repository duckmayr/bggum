#include "../inst/include/ggumR.h"

using namespace Rcpp;

//' GGUM getPrior
//' 
//' Get the prior probability of values for the GGUM parameters theta, alpha,
//' delta, and tau.
//' 
//' @param cv A numeric vector of length one; the current value of the
//'   parameter of interest
//' 
//' @return A numeric vector of length one. The prior probability of observing
//'   \code{cv} for the paramenter of interest.
//' @rdname getPrior
//' @export
//[[Rcpp::export]]
double getPriorTheta(double cv){
    return R::dnorm(cv, 0, 1, 0);
}

//' @export
//[[Rcpp::export]]
double getPriorAlpha(double cv){
    return d_4beta(cv, 1.5, 1.5, 0.25, 4);
}

//' @export
//[[Rcpp::export]]
double getPriorDelta(double cv){
    return d_4beta(cv, 2, 2, -5, 5);
}

//' @export
//[[Rcpp::export]]
double getPriorTaus(double cv){
    return d_4beta(cv, 2, 2, -6, 6);
}
