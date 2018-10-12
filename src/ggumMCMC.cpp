#include "bggum.h"

using namespace Rcpp;

//[[Rcpp::export(.ggumMCMC)]]
NumericMatrix ggumMCMC(IntegerMatrix data, int n, int m, int iterations,
        int burn_iterations, int flip_interval,
        NumericVector thetas, NumericVector alphas,
        NumericVector deltas, List taus, IntegerVector K,
        double th_prior_mean, double th_prior_sd, double a_shape1,
        double a_shape2, double a_a, double a_b, double d_shape1,
        double d_shape2, double d_a, double d_b, double t_shape1,
        double t_shape2, double t_a, double t_b, List SDs){
    NumericVector theta_SDs = as<NumericVector>(SDs[0]);
    NumericVector alpha_SDs = as<NumericVector>(SDs[1]);
    NumericVector delta_SDs = as<NumericVector>(SDs[2]);
    NumericVector tau_SDs = as<NumericVector>(SDs[3]);
    // First run the burn in
    // Set up the progress display
    Rcout.precision(1);
    double adv_prog = 10000.0 / burn_iterations;
    int current_break = 1;
    if ( burn_iterations > 0 ) {
        Rcout << "\rBurning in:          0%";
    }
    // Then we run the burn in iterations:
    for ( int iter = 0; iter < burn_iterations; ++iter ) {
        if ( (iter+1) % 100 == 0 ) {
            // Every 100 iterations we check for a user interrupt
            // and update the progress bar
            checkUserInterrupt();
            double current_prog = current_break * adv_prog;
            Rcout << "\rBurning in:          " << std::fixed << current_prog << "%";
            current_break += 1;
        }
        // Then we update the variables
        for ( int i = 0; i < n; ++i ) {
            // Copy the parameter of interest:
            double theta = thetas[i];
            // Replace it (or not):
            thetas[i] = update_theta_MCMC(data(i, _), theta, alphas,
                    deltas, taus, theta_SDs[i], th_prior_mean, th_prior_sd);
        }
        for ( int j = 0; j < m; ++j ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            alphas[j] = update_alpha_MCMC(data(_, j), thetas, alpha,
                    delta, taus[j], alpha_SDs[j], a_shape1, a_shape2, a_a, a_b);
        }
        for ( int j = 0; j < m; ++j ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            deltas[j] = update_delta_MCMC(data(_, j), thetas, alpha,
                    delta, taus[j], delta_SDs[j], d_shape1, d_shape2, d_a, d_b);
        }
        for ( int j = 0; j < m; ++j ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            NumericVector thisTau = taus[j];
            for ( int k = 1; k < K[j]; ++k ) {
                thisTau[k] = update_tau_MCMC(k, data(_, j), thetas,
                        alpha, delta, thisTau, tau_SDs[j], t_shape1, t_shape2,
                        t_a, t_b);
            }
            taus[j] = thisTau;
        }
    }
    if ( burn_iterations > 0 ) {
        Rcout << "\rBurning in:          100.0%\n";
    }
    // This makes an empty matrix to store parameter values for every iteration:
    NumericMatrix chainMatrix(iterations, n+(2*m)+sum(K-1));
    // set up progress display
    adv_prog = 10000.0 / iterations;
    current_break = 1;
    Rcout << "\rSampling posterior:  0%";
    // Then we run the MCMC sampler:
    for ( int iter = 0; iter < iterations; ++iter ) {
        if ( (iter+1) % 100 == 0 ) {
            // Every 100 iterations we check for a user interrupt
            // and update the progress bar
            checkUserInterrupt();
            double current_prog = current_break * adv_prog;
            Rcout << "\rSampling posterior:  " << std::fixed << current_prog << "%";
            current_break += 1;
        }
        if ( (iter+1) % flip_interval == 0 ) {
            thetas = thetas * -1.0;
            deltas = deltas * -1.0;
        }
        // Then we update the variables
        for ( int i = 0; i < n; ++i ) {
            // Copy the parameter of interest:
            double theta = thetas[i];
            // Replace it (or not):
            thetas[i] = update_theta_MCMC(data(i, _), theta, alphas,
                    deltas, taus, theta_SDs[i], th_prior_mean, th_prior_sd);
            // And store the parameter value for this iteration in the matrix:
            chainMatrix(iter, i) = thetas[i];
        }
        for ( int j = 0; j < m; ++j ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            alphas[j] = update_alpha_MCMC(data(_, j), thetas, alpha,
                    delta, taus[j], alpha_SDs[j], a_shape1, a_shape2, a_a, a_b);
            chainMatrix(iter, j+n) = alphas[j];
        }
        for ( int j = 0; j < m; ++j ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            deltas[j] = update_delta_MCMC(data(_, j), thetas, alpha,
                    delta, taus[j], delta_SDs[j], d_shape1, d_shape2, d_a, d_b);
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
                thisTau[k] = update_tau_MCMC(k, data(_, j), thetas,
                        alpha, delta, thisTau, tau_SDs[j], t_shape1, t_shape2,
                        t_a, t_b);
                chainMatrix(iter, n+(2*m)+Ksum+k-1) = thisTau[k];
            }
            taus[j] = thisTau;
            Ksum += K[j]-1;
        }
    }
    // Close out the progress display
    Rcout << "\rSampling posterior:  100.0%\n";
    // And finally, we return the chain matrix
    return chainMatrix;
}
