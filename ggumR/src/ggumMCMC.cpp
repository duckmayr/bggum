#include "../inst/include/ggumR.h"

using namespace Rcpp;

//' GGUM MCMC Sampler
//'
//' MCMC sampler for the generalized graded unfolding model (GGUM), utilizing
//' a Metropolis-Hastings algorithm
//'
//' @param responseMatrix A numeric matrix giving the response by each
//'   respondent to each item
//' @param Kvector A numeric vector of length m (the number of items), each
//'   element j of which gives the number of options (K) for item j
//' @param iterations A numeric vector of length one; the number of iterations
//'
//' @return A chain matrix; a numeric matrix with \code{iterations} rows and
//'   one column for every parameter of the model, so that each element of the
//'   matrix gives the value of a parameter for a particular iteration of the
//'   algorithm
//' @export
//[[Rcpp::export]]
NumericMatrix ggumMCMC(NumericMatrix responseMatrix, IntegerVector Kvector,
                       int iterations){
    int n = responseMatrix.nrow();
    int m = responseMatrix.ncol();
    NumericVector thetas = rnorm(n, 0, 1);
    NumericVector alphas = r4beta(m, 1.5, 1.5, 0.25, 4);
    NumericVector deltas = r4beta(m,  2, 2, -5, 5);
    List taus(m);
    for ( int j = 0; j < m; j++ ){
        int K = Kvector[j];
        NumericVector thisTau(K);
        thisTau[Range(1, K-1)] = r4beta(K-1, 2, 2, -6, 6);
        taus[j] = thisTau;
    }
    int numberOfTaus = 0;
    for ( int j = 0; j < m; j++ ) {
        numberOfTaus += Kvector[j];
    }
    NumericMatrix chainMatrix(iterations, n+(2*m)+numberOfTaus);
    for ( int iter = 0; iter < 1000; iter++ ) {
        for ( int i = 0; i < n; i++ ) {
            double theta = thetas[i];
            thetas[i] = acceptanceTheta(responseMatrix(i, _), theta, alphas,
                    deltas, taus, 1);
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
    for ( int iter = 1000; iter < iterations; iter ++ ) {
        for ( int i = 0; i < n; i++ ) {
            double theta = thetas[i];
            NumericVector pastDraws = chainMatrix(_, i);
            double SD = sd(pastDraws[Range(iter-1000, iter)]);
            thetas[i] = acceptanceTheta(responseMatrix(i, _), theta, alphas,
                    deltas, taus, SD);
            chainMatrix(iter, i) = thetas[i];
        }
        for ( int j = 0; j < m; j++ ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            NumericVector pastDraws = chainMatrix(_, j+n);
            double SD = sd(pastDraws[Range(iter-1000, iter)]);
            alphas[j] = acceptanceAlpha(responseMatrix(_, j), thetas, alpha,
                    delta, taus[j], SD);
            chainMatrix(iter, j+n) = alphas[j];
        }
        for ( int j = 0; j < m; j++ ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            NumericVector pastDraws = chainMatrix(_, j+n+m);
            double SD = sd(pastDraws[Range(iter-1000, iter)]);
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
                double SD = sd(pastDraws[Range(iter-1000, iter)]);
                thisTau[k] = acceptanceTau(k, responseMatrix(_, j), thetas,
                        alpha, delta, thisTau, SD);
                chainMatrix(iter, n+(2*m)+Ksum+k) = thisTau[k];
            }
            taus[j] = thisTau;
            Ksum += K;
        }
    }
    return chainMatrix;
}
