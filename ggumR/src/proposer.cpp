#include "../inst/include/ggumR.h"

using namespace Rcpp;

//' GGUM Proposer
//' 
//' Propose new values for the GGUM parameters Theta, Alpha, Delta, and Tau
//' 
//' Given a current value of the parameter of interest,
//' and a sigma parameter for the proposal distribution,
//' a new proposal for the parameter's value is given.
//' The following proposal densities are used:
//' \itemize{
//'   \item Theta -- A location-scale T distribution shifted by \code{cv} and
//'     scaled by \code{SD}, truncated at -10 and 10
//'   \item Alpha -- A location-scale T distribution shifted by \code{cv} and
//'     scaled by \code{SD}, truncated at 0.25 and 4
//'   \item Delta -- A location-scale T distribution shifted by \code{cv} and
//'     scaled by \code{SD}, truncated at -5 and 5
//'   \item Tau -- A location-scale T distribution shifted by \code{cv} and
//'     scaled by \code{SD}, truncated at -6 and 6
//' }
//'
//' @param cv A numeric vector of length one;
//'   the current value of the parameter of interest.
//' @param SD A numeric vector of length one;
//'   the sigma parameter used to generate the proposed value.
//' 
//' @return A numeric vector of length one;
//'   the proposed value for the parameter of interest.
//' @name GGUM MCMC Proposal Densities
//' @aliases proposerTheta, proposerAlpha, proposerDelta, proposerTau,
//'   proposer
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
