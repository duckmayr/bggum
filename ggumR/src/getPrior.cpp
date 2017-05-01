#include "../inst/include/ggumR.h"

using namespace Rcpp;

//' GGUM getPrior
//' 
//' Get the prior probability of values for the GGUM parameters theta, alpha,
//' delta, and tau.
//'
//' Following de la Torre et al (2006), we use the following prior
//' distributions for MCMC estimation of GGUM parameters:
//' \itemize{
//'   \item Theta -- The standard normal
//'   \item Alpha -- A four parameter beta distribution, with shape parameters
//'     1.5, 1.5, minimum value 0.25, and maximum value 4
//'   \item Delta -- A four parameter beta distribution, with shape parameters
//'     2, 2, minimum value -5, and maximum value 5
//'   \item Tau -- A four parameter beta distribution, with shape parameters
//'     2, 2, minimum value -6, and maximum value 6
//' }
//' 
//' @param cv A numeric vector of length one; the current value of the
//'   parameter of interest
//' 
//' @return A numeric vector of length one. The prior probability of observing
//'   \code{cv} for the paramenter of interest.
//'
//' @references de la Torre, Jimmy, Stephen Stark, and Oleksandr S.
//'   Chernyshenko. 2006. ``Markov Chain Monte Carlo Estimation of Item
//'   Parameters for the Generalized Graded Unfolding Model." \emph{Applied
//'   Psychological Measurement} 30(3): 216--232.
//'
//' @name GGUM Priors
//' @aliases getPriorTheta, getPriorAlpha, getPriorDelta, getPriorTau, getPrior
//' @rdname getPrior
//' @export
//[[Rcpp::export]]
double getPriorTheta(double cv){
    return R::dnorm(cv, 0, 1, 0);
}

//' @rdname getPrior
//' @export
//[[Rcpp::export]]
double getPriorAlpha(double cv){
    return d_4beta(cv, 1.5, 1.5, 0.25, 4);
}

//' @rdname getPrior
//' @export
//[[Rcpp::export]]
double getPriorDelta(double cv){
    return d_4beta(cv, 2, 2, -5, 5);
}

//' @rdname getPrior
//' @export
//[[Rcpp::export]]
double getPriorTaus(double cv){
    return d_4beta(cv, 2, 2, -6, 6);
}
