#include "../inst/include/ggum.h"

using namespace Rcpp;

//[[Rcpp::export(.ggumMCMC)]]
NumericMatrix ggumMCMC(IntegerMatrix responseMatrix, int iterations,
        int burn_iterations, NumericVector thetas, NumericVector alphas,
        NumericVector deltas, List taus, IntegerVector K, int n, int m){
    // First we run the burn-in
    // (for now, this will automatically tune the proposals)
    List SDs = tune_proposals(responseMatrix, thetas, alphas, deltas, taus, K,
                              burn_iterations, n, m);
    NumericVector theta_SDs = as<NumericVector>(SDs[0]);
    NumericVector alpha_SDs = as<NumericVector>(SDs[1]);
    NumericVector delta_SDs = as<NumericVector>(SDs[2]);
    NumericVector tau_SDs = as<NumericVector>(SDs[3]);
    // This makes an empty matrix to store parameter values for every iteration:
    NumericMatrix chainMatrix(iterations, n+(2*m)+sum(K));
    // set up progress display
    Rcout.precision(1);
    double adv_prog = 10000.0 / iterations;
    int current_break = 1;
    Rcout << "\rRunning sampler: 0%";
    // Then we run the MCMC sampler:
    for ( int iter = 0; iter < iterations; ++iter ) {
        if ( (iter+1) % 100 == 0 ) {
            // Every 100 iterations we check for a user interrupt
            // and update the progress bar
            checkUserInterrupt();
            double current_prog = current_break * adv_prog;
            Rcout << "\rRunning sampler: " << std::fixed << current_prog << "%";
            current_break += 1;
        }
        // Then we update the variables
        for ( int i = 0; i < n; ++i ) {
            // Copy the parameter of interest:
            double theta = thetas[i];
            // Replace it (or not):
            thetas[i] = update_theta_MCMC(responseMatrix(i, _), theta, alphas,
                    deltas, taus, theta_SDs[i]);
            // And store the parameter value for this iteration in the matrix:
            chainMatrix(iter, i) = thetas[i];
        }
        for ( int j = 0; j < m; ++j ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            alphas[j] = update_alpha_MCMC(responseMatrix(_, j), thetas, alpha,
                    delta, taus[j], alpha_SDs[j]);
            chainMatrix(iter, j+n) = alphas[j];
        }
        for ( int j = 0; j < m; ++j ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            deltas[j] = update_delta_MCMC(responseMatrix(_, j), thetas, alpha,
                    delta, taus[j], delta_SDs[j]);
            chainMatrix(iter, j+n+m) = deltas[j];
        }
        // For taus, we need to keep track of where in the chain matrix we
        // should be storing values, so we create the variable Ksum, which
        // when added with (n + 2 * m + k) gives us the column of the matrix
        // in which to inser values
        int Ksum = 0;
        for ( int j = 0; j < m; ++j ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            NumericVector thisTau = taus[j];
            for ( int k = 1; k < K[j]; ++k ) {
                thisTau[k] = update_tau_MCMC(k, responseMatrix(_, j), thetas,
                        alpha, delta, thisTau, tau_SDs[j]);
                chainMatrix(iter, n+(2*m)+Ksum+k) = thisTau[k];
            }
            taus[j] = thisTau;
            Ksum += K[j];
        }
    }
    // Close out the progress display
    Rcout << "\n";
    // And finally, we return the chain matrix
    return chainMatrix;
}
