#include "../inst/include/ggumR.h"

using namespace Rcpp;

//' GGUM Proposer
//' 
//' Propose new values for the parameters Theta, Alpha, Delta, and Taus
//' for the GGUM.
//' 
//' @param cv A numeric vector of length one. The current value of the parameter of interest.
//' @param SD A numeric vector of length one. The standard deviation used to generate the 
//' proposed value.
//' 
//' @return A numeric vector of length one. The propose value for the paramenter of interest.
//' @rdname ggumProposer
//' @export
//[[Rcpp::export]]
double proposerTheta(NumericVector cv, NumericVector SD){
   double pv = r_trunclst(1, cv[0], SD[0], -10, 10);
   return pv;
}

//' @rdname ggumProposer
//' @export
//[[Rcpp::export]]
double proposerAlpha(NumericVector cv, NumericVector SD) {
   double pv = r_trunclst(1, cv[0], SD[0], 0.25, 4);
   return pv;
}

//' @rdname ggumProposer
//' @export
//[[Rcpp::export]]
double proposerDelta(NumericVector cv, NumericVector SD){
   double pv = r_trunclst(1, cv[0], SD[0], -5, 5);
   return pv;
}

//' @rdname ggumProposer
//' @export
//[[Rcpp::export]]
double proposerTau(NumericVector cv, NumericVector SD){
   double pv = r_trunclst(1, cv[0], SD[0], -6, 6);
   return pv;
}