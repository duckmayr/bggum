#include "model.h" 

// [[Rcpp::export(.ggumMCMC)]]
Rcpp::NumericMatrix ggumMCMC(Rcpp::IntegerMatrix data_,
                             Rcpp::NumericMatrix theta_,
                             Rcpp::NumericMatrix alpha_,
                             Rcpp::NumericMatrix delta_,
                             Rcpp::List tau_,
                             Rcpp::IntegerVector K_,
                             Rcpp::List SDs,
                             Rcpp::NumericVector temps_,
                             Rcpp::NumericVector alpha_parameters_,
                             Rcpp::NumericVector delta_parameters_,
                             Rcpp::NumericVector tau_parameters_,
                             int sample_iterations,
                             int burn_iterations,
                             int flip_interval) {
    // Separate the proposal s.d. list
    Rcpp::NumericVector theta_SDs = get_vec_from_list(SDs, 0);
    Rcpp::NumericVector alpha_SDs = get_vec_from_list(SDs, 1);
    Rcpp::NumericVector delta_SDs = get_vec_from_list(SDs, 2);
    Rcpp::NumericVector tau_SDs   = get_vec_from_list(SDs, 3);
    // Set up the Model object
    Model mod(data_, theta_, alpha_, delta_, tau_, K_,
              theta_SDs, alpha_SDs, delta_SDs, tau_SDs,
              alpha_parameters_, delta_parameters_, tau_parameters_,
              temps_);
    // Set up a matrix to store the posterior draws
    int n_params = mod.n + 2 * mod.m + Rcpp::sum(mod.K - 1);
    Rcpp::NumericMatrix result(sample_iterations, n_params);
    // Set up the progress display
    double progress_increment = (10000.0 / burn_iterations);
    double progress = 0.0;
    if ( burn_iterations > 0 ) {
        Rprintf("Burning in:          %7.3f %%", progress);
    }
    // Then we run the burn in iterations:
    for ( int iter = 0; iter < burn_iterations; ++iter ) { // for each iteration
        // Update parameters
        mod.update_theta(0);
        mod.update_alpha(0);
        mod.update_delta(0);
        mod.update_tau(0);
        if ( (iter+1) % 100 == 0 ) {
            // Every 100 iterations, check for user interrupt & update progress
            Rcpp::checkUserInterrupt();
            progress += progress_increment;
            Rprintf("\rBurning in:          %7.3f %%", progress);
        }
    }
    if ( burn_iterations > 0 ) {
        Rprintf("\rBurning in:          %7.3f %%\n", progress);
    }
    // set up progress display
    progress_increment = (10000.0 / sample_iterations);
    progress = 0.0;
    Rprintf("Sampling posterior:  %7.3f %%", progress);
    // now run the sampler
    for ( int iter = 0; iter < sample_iterations; ++iter ) {
        // Update parameters
        mod.update_theta(0);
        mod.update_alpha(0);
        mod.update_delta(0);
        mod.update_tau(0);
        if ( (iter+1) % 100 == 0 ) {
            Rcpp::checkUserInterrupt();
            progress += progress_increment;
            Rprintf("\rSampling posterior:  %7.3f %%", progress);
        }
        if ( (iter+1) % flip_interval == 0 ) {
            mod.theta = mod.theta * -1.0;
            mod.delta = mod.delta * -1.0;
        }
        for ( int i = 0; i < mod.n; ++i ) {
            result(iter, i) = mod.theta(i, 0);
        }
        int K_ind = 0;
        int idx = 0;
        for ( int j = 0; j < mod.m; ++j ) {
            idx = mod.n + j;
            result(iter, idx) = mod.alpha(j, 0);
            idx += mod.m;
            result(iter, idx) = mod.delta(j, 0);
            idx += (mod.m + K_ind - 1);
            for ( int k = 1; k < mod.K[j]; ++k ) {
                idx += 1;
                result(iter, idx) = mod.tau[0][j][k];
            }
            K_ind += (mod.K[j] - 2);
        }
    }
    // return the cold chain
    Rprintf("\rSampling posterior:  %7.3f %%\n", progress);
    return result;
}

