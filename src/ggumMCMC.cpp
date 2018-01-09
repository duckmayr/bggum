#include "../inst/include/ggum.h"

using namespace Rcpp;

//' GGUM MCMC Sampler
//'
//' MCMC sampler for the generalized graded unfolding model (GGUM), utilizing
//' a Metropolis-Hastings algorithm
//'
//' \code{ggumMCMC} provides \code{R} implementation of an MCMC sampler for
//' the GGUM, based heavily on the algorithm given in de la Torre et al (2006);
//' though the package allows parameter estimation from \code{R},
//' the functions are actually written in \code{C++} to allow for reasonable
//' execution time.
//' Our sampler creates random ititial values for the parameters of the model,
//' according to their prior distributions (see \code{\link{getPrior}}).
//' Then, for the first 5000 iterations, the sampler, one parameter at a time,
//' will make a proposal from the truncated location-scale T distribution
//' specified in \code{\link{proposer}} where \eqn{SD = 1}, and accept the
//' proposal probabilistically according to the ratio in
//' \code{\link{acceptance}}. For the remainder of the iterations, the same
//' process is followed, but the \eqn{\sigma} parameter for the proposal
//' density is the standard deviation of the previous 5000 values of the
//' parameter. A matrix is returned giving the value of every parameter at
//' every iteration.
//'
//' @section Warning:
//' Typically, up to 5000 iterations are needed for convergence,
//' and specifying less than 5000 iterations will cause an error.
//' Also, theta and delta parameters often have multi-modal posterior
//' distributions, and in some cases chains may converge to the wrong
//' mode for both parameters; in other words, the questions with the most
//' negative delta parameters will be estimated as having the most positive
//' delta parameters and vice versa, but only if the same happens with
//' the theta parameters. Future package updates will provide a way to
//' specify the expected sign of at least one theta and delta parameter
//' to prevent this occurrence.
//'
//' @param responseMatrix A numeric matrix giving the response by each
//'   respondent to each item
//' @param Kvector A numeric vector of length m (the number of items), each
//'   element j of which gives the number of options (K) for item j
//' @param iterations A numeric vector of length one; the number of iterations
//'   (NOTE: \code{iterations} should be at least 10000, and preferably around
//'   25000, though only values of 5000 or less will cause an error)
//' @param low A numeric vector of length one giving the row number for a
//'   respondent whose theta parameter will be restricted to be negative
//' @param high A numeric vector of length one giving the row number for a
//'   respondent whose theta parameter will be restricted to be positive
//'
//' @return A chain matrix; a numeric matrix with \code{iterations} rows and
//'   one column for every parameter of the model, so that each element of the
//'   matrix gives the value of a parameter for a particular iteration of the
//'   MCMC algorithm.
//'
//' @seealso \code{\link{ggumProbability}}, \code{\link{acceptance}},
//'   \code{\link{getPrior}}, \code{\link{proposer}}
//'
//' @references Roberts, James S., John R. Donoghue, and James E. Laughlin.
//'   2000. ``A General Item Response Theory Model for Unfolding
//'   Unidimensional Polytomous Responses." \emph{Applied Psychological
//'   Measurement} 24(1): 3--32.
//' @references de la Torre, Jimmy, Stephen Stark, and Oleksandr S.
//'   Chernyshenko. 2006. ``Markov Chain Monte Carlo Estimation of Item
//'   Parameters for the Generalized Graded Unfolding Model." \emph{Applied
//'   Psychological Measurement} 30(3): 216--232.
//'   algorithm
//' @export
//[[Rcpp::export]]
NumericMatrix ggumMCMC(IntegerMatrix responseMatrix, IntegerVector Kvector,
                       int iterations, int low, int high){
    // This deals with indexing differences between R and C++
    low -= 1;
    high -= 1;
    // It will be convenient to store the number of items and respondents:
    int n = responseMatrix.nrow();
    int m = responseMatrix.ncol();
    // Then we draw initial parameter values from their prior distributions:
    NumericVector thetas = rnorm(n, 0, 1);
    NumericVector alphas = r4beta(m, 1.5, 1.5, 0.25, 4);
    NumericVector deltas = r4beta(m,  2, 2, -5, 5);
    List taus(m);
    for ( int j = 0; j < m; j++ ){
        int K = Kvector[j];
        NumericVector thisTau(K);
        thisTau[Range(1, K-1)] = r4beta(K-1, 2, 2, -2, 0);
        taus[j] = thisTau;
    }
    // It will also be convenient to store the total number of tau parameters:
    int numberOfTaus = 0;
    for ( int j = 0; j < m; j++ ) {
        numberOfTaus += Kvector[j];
    }
    // This makes an empty matrix to store parameter values for every iteration:
    NumericMatrix chainMatrix(iterations, n+(2*m)+numberOfTaus);
    // Now we need to deal with restrictions;
    // if our theta draws for don't respect the restrictions,
    // we need to change them:
    if ( thetas[low] > 0 ) {
        thetas[low] = thetas[low] * -1;
    }
    if ( thetas[high] < 0 ) {
        thetas[high] = thetas[high] * -1;
    }
    // Now it will be convenient to have a vector giving the indices
    // of the unconstrained theta parameters (there are n - 2 of them):
    IntegerVector unconstrained(n-2);
    // If one of the constrained parameters is theta[0],
    // we need to make sure the first element of this vector is 1:
    if ( low == 0 || high == 0 ) {
        unconstrained[0] = 1;
    }
    // Then we can get the rest of the vector's elements easily:
    for ( int i = 1; i < n-2; i++ ) {
        unconstrained[i] = unconstrained[i-1] + 1;
        if ( unconstrained[i] == low || unconstrained[i] == high ) {
            unconstrained[i] += 1;
        }
    }
    // For the first 5000 iterations we use a fixed sigma for proposals:
    for ( int iter = 0; iter < 5000; iter++ ) {
        double theta = thetas[low];
        thetas[low] = acceptanceThetaNeg(responseMatrix(low, _), theta,
                alphas, deltas, taus, 1);
        chainMatrix(iter, low) = thetas[low];
        theta = thetas[high];
        thetas[high] = acceptanceThetaPos(responseMatrix(high, _), theta,
                alphas, deltas, taus, 1);
        chainMatrix(iter, high) = thetas[high];
        for ( int ind = 0; ind < n-2; ind++ ) {
            int i = unconstrained[ind];
            // Copy the parameter of interest:
            double theta = thetas[i];
            // Replace it (or not):
            thetas[i] = acceptanceTheta(responseMatrix(i, _), theta, alphas,
                    deltas, taus, 1);
            // And store the parameter value for this iteration in the matrix:
            chainMatrix(iter, i) = thetas[i];
        }
        for ( int j = 0; j < m; j++ ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            alphas[j] = acceptanceAlpha(responseMatrix(_, j), thetas, alpha,
                    delta, taus[j], 1);
            chainMatrix(iter, j+n) = alphas[j];
        }
        for ( int j = 0; j < m; j++ ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            deltas[j] = acceptanceDelta(responseMatrix(_, j), thetas, alpha,
                    delta, taus[j], 1);
            chainMatrix(iter, j+n+m) = deltas[j];
        }
        // For taus, we need to keep track of where in the chain matrix we
        // should be storing values, so we create the variable Ksum, which
        // when added with (n + 2 * m + k) gives us the column of the matrix
        // in which to inser values
        int Ksum = 0;
        for ( int j = 0; j < m; j++ ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            int K = Kvector[j];
            NumericVector thisTau = taus[j];
            for ( int k = 1; k < K; k++ ) {
                thisTau[k] = acceptanceTau(k, responseMatrix(_, j), thetas,
                        alpha, delta, thisTau, 1);
                chainMatrix(iter, n+(2*m)+Ksum+k) = thisTau[k];
            }
            taus[j] = thisTau;
            Ksum += K;
        }
    }
    // For the remaining iterations we use the standard deviation of the past
    // 5000 draws for the sigma parameter for every parameter
    for ( int iter = 5000; iter < iterations; iter ++ ) {
        if ( iter % 1000 == 0 ) {
            checkUserInterrupt();
        }
        double theta = thetas[low];
        NumericVector pastDraws = chainMatrix(_, low);
        double SD = sd(pastDraws[Range(iter-5000, iter)]);
        thetas[low] = acceptanceThetaNeg(responseMatrix(low, _), theta,
                alphas, deltas, taus, SD);
        chainMatrix(iter, low) = thetas[low];
        theta = thetas[high];
        pastDraws = chainMatrix(_, high);
        SD = sd(pastDraws[Range(iter-5000, iter)]);
        thetas[high] = acceptanceThetaPos(responseMatrix(high, _), theta,
                alphas, deltas, taus, 1);
        chainMatrix(iter, high) = thetas[high];
        for ( int ind = 0; ind < n-2; ind++ ) {
            int i = unconstrained[ind];
            double theta = thetas[i];
            NumericVector pastDraws = chainMatrix(_, i);
            double SD = sd(pastDraws[Range(iter-5000, iter)]);
            thetas[i] = acceptanceTheta(responseMatrix(i, _), theta, alphas,
                    deltas, taus, SD);
            chainMatrix(iter, i) = thetas[i];
        }
        for ( int j = 0; j < m; j++ ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            NumericVector pastDraws = chainMatrix(_, j+n);
            double SD = sd(pastDraws[Range(iter-5000, iter)]);
            alphas[j] = acceptanceAlpha(responseMatrix(_, j), thetas, alpha,
                    delta, taus[j], SD);
            chainMatrix(iter, j+n) = alphas[j];
        }
        for ( int j = 0; j < m; j++ ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            NumericVector pastDraws = chainMatrix(_, j+n+m);
            double SD = sd(pastDraws[Range(iter-5000, iter)]);
            deltas[j] = acceptanceDelta(responseMatrix(_, j), thetas, alpha,
                    delta, taus[j], SD);
            chainMatrix(iter, j+n+m) = deltas[j];
        }
        int Ksum = 0;
        for ( int j = 0; j < m; j++ ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            int K = Kvector[j];
            NumericVector thisTau = taus[j];
            for ( int k = 1; k < K; k++ ) {
                NumericVector pastDraws = chainMatrix(_, n+Ksum+k+(2*m));
                double SD = sd(pastDraws[Range(iter-5000, iter)]);
                thisTau[k] = acceptanceTau(k, responseMatrix(_, j), thetas,
                        alpha, delta, thisTau, SD);
                chainMatrix(iter, n+(2*m)+Ksum+k) = thisTau[k];
            }
            taus[j] = thisTau;
            Ksum += K;
        }
    }
    // And finally, we return the chain matrix
    return chainMatrix;
}
