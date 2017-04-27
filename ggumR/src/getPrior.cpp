#include "../inst/include/ggumR.h"

using namespace Rcpp;

//' GGUM getPrior
//' 
//' Get the prior values for the parameters Theta, Alpha, Delta, and Taus
//' for the GGUM.
//' 
//' @param cv A numeric vector of length one. The current value of the parameter of interest
//' @param tau The vector of current values of the items' tau parameters
//' 
//' @return A numeric vector of length one. The propose value for the paramenter of interest.
//' @rdname getPrior
//' @export
//[[Rcpp::export]]
double getPriorTheta(double cv){
  double cvPrior = d_4beta(cv, 2, 2, -5, 5);
}

//' @export
//[[Rcpp::export]]
double getPriorAlpha(double cv){
  double cvPrior = d_4beta(cv, 2, 2, -5, 5);
}

//' @export
//[[Rcpp::export]]
double getPriorDelta(double cv){
  double cvPrior = d_4beta(cv, 2, 2, -5, 5);
}

//' @export
//[[Rcpp::export]]
double getPriorTaus(double cv, NumericVector taus){
  double cvPrior = d_4beta(cv, 2, 2, -5, 5);
  double cvPrior = d_4beta(taus[k], 2, 2, -6, 6);
}