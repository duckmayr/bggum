#include "model.h" 

// [[Rcpp::export(.tune_temperatures)]]
Rcpp::NumericVector tune_temps(Rcpp::IntegerMatrix data_,
                               Rcpp::NumericMatrix theta_,
                               Rcpp::NumericMatrix alpha_,
                               Rcpp::NumericMatrix delta_,
                               Rcpp::List tau_,
                               Rcpp::IntegerVector K_,
                               Rcpp::List SDs,
                               Rcpp::NumericVector alpha_parameters_,
                               Rcpp::NumericVector delta_parameters_,
                               Rcpp::NumericVector tau_parameters_,
                               int n_temps,
                               int iterations,
                               int n_draws) {
    // Separate the proposal s.d. list
    Rcpp::NumericVector theta_SDs = get_vec_from_list(SDs, 0);
    Rcpp::NumericVector alpha_SDs = get_vec_from_list(SDs, 1);
    Rcpp::NumericVector delta_SDs = get_vec_from_list(SDs, 2);
    Rcpp::NumericVector tau_SDs   = get_vec_from_list(SDs, 3);
    // Set up the Model object
    Rcpp::NumericVector temps(n_temps, 1.0);
    Model mod(data_, theta_, alpha_, delta_, tau_, K_,
              theta_SDs, alpha_SDs, delta_SDs, tau_SDs,
              alpha_parameters_, delta_parameters_, tau_parameters_,
              temps);
    // Set up storage for temperature values at each iteration
    Rcpp::NumericVector tmp(n_draws);
    // Set up algorithm variables
    double rho1 = -2.0;
    double alph = 0.0;
    double z = 1.0; // what Atchade et al. call n
    mod.temps[1] = 1.0 / (1.0 + exp(rho1));
    // Set up progress display
    double progress_increment = 100.0 / ((iterations / 100.0) * (n_temps - 1));
    double progress = 0.0;
    Rprintf("Tuning temperatures: %7.3f %%", progress);
    // Run algorithm
    for ( int t = 1; t < n_temps; ++t ) {
        for ( int iter = 0; iter < iterations; ++iter ) {
            if ( (iter+1) % 100 == 0 ) {
                // Every 100 iterations we check for a user interrupt
                // and update the progress bar
                Rcpp::checkUserInterrupt();
                progress += progress_increment;
                Rprintf("\rTuning temperatures: %7.3f %%", progress);
            }
            for ( int l = 0; l < 2; ++l ) {
                mod.update_theta(t - l);
                mod.update_alpha(t - l);
                mod.update_delta(t - l);
                mod.update_tau(t - l);
            }
            alph = mod.temps[t - 1] - mod.temps[t];
            double chain1_prior = 0.0, chain2_prior = 0.0;
            for ( int i = 0; i < mod.n; ++i ) {
                chain1_prior += mod.theta_prior(mod.theta(i, t - 1));
                chain2_prior += mod.theta_prior(mod.theta(i, t));
            }
            for ( int j = 0; j < mod.m; ++j ) {
                chain1_prior += mod.alpha_prior(mod.alpha(j, t - 1));
                chain1_prior += mod.delta_prior(mod.delta(j, t - 1));
                chain2_prior += mod.alpha_prior(mod.alpha(j, t));
                chain2_prior += mod.delta_prior(mod.delta(j, t));
                for ( int k = 1; k < mod.K[j]; ++k ) {
                    chain1_prior += mod.tau_prior(mod.tau[t - 1][j][k]);
                    chain2_prior += mod.tau_prior(mod.tau[t][j][k]);
                }
            }
            double chain1_log_likelihood = mod.log_likelihood(t - 1);
            double chain2_log_likelihood = mod.log_likelihood(t);
            double multiplier = chain2_log_likelihood - chain1_log_likelihood;
            multiplier += (chain2_prior - chain1_prior);
            alph *= multiplier;
            alph = exp(alph);
            if ( alph > 1.0 ) {
                alph = 1.0;
            }
            rho1 += (1.0 / z) * (alph - 0.2338071);
            mod.temps[t] = mod.temps[t - 1] / (1.0 + exp(rho1));
            if ( iter >= (iterations - n_draws) ) {
                tmp[iterations - iter - 1] = mod.temps[t];
            }
        }
        z += 1.0;
        rho1 = -2.0;
        mod.temps[t] = Rcpp::mean(tmp);
        if ( t < (n_temps-1) ) {
            mod.temps[t + 1] = mod.temps[t] / (1.0 + exp(rho1));
        }
    }
    Rprintf("\rTuning temperatures: %7.3f %%\n", progress);
    return mod.temps;
}

