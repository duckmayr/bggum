#include "model.h"

// [[Rcpp::export(.tune_proposals)]]
Rcpp::List tune_proposals(Rcpp::IntegerMatrix data_,
                          Rcpp::NumericMatrix theta_,
                          Rcpp::NumericMatrix alpha_,
                          Rcpp::NumericMatrix delta_,
                          Rcpp::List tau_,
                          Rcpp::IntegerVector K_,
                          Rcpp::NumericVector alpha_parameters_,
                          Rcpp::NumericVector delta_parameters_,
                          Rcpp::NumericVector tau_parameters_,
                          Rcpp::NumericVector temps_,
                          int tune_iterations) {
    // set up progress display
    double progress_increment = 1000.0 / tune_iterations;
    double progress = 0.0;
    Rprintf("Tuning proposals:    %7.3f %%", progress);
    // set up vectors for proposal SDs
    Rcpp::NumericVector theta_SD(data_.nrow(), 0.01);
    Rcpp::NumericVector alpha_SD(data_.ncol(), 0.01);
    Rcpp::NumericVector delta_SD(data_.ncol(), 0.01);
    Rcpp::NumericVector tau_SD(data_.ncol(), 0.01);
    // Set up the Model object
    Model mod(data_, theta_, alpha_, delta_, tau_, K_,
              theta_SD, alpha_SD, delta_SD, tau_SD,
              alpha_parameters_, delta_parameters_, tau_parameters_,
              temps_);
    // set up vectors to store number of acceptances per cycle
    Rcpp::IntegerVector theta_accepts(mod.n);
    Rcpp::IntegerVector alpha_accepts(mod.m);
    Rcpp::IntegerVector delta_accepts(mod.m);
    Rcpp::IntegerVector tau_accepts(mod.m);
    for ( int iter = 0; iter < tune_iterations; ++iter ) {
        if ( (iter+1) % 10 == 0 ) {
            // Check for user interruption
            Rcpp::checkUserInterrupt();
            // Update progress display
            progress += progress_increment;
            Rprintf("\rTuning proposals:    %7.3f %%", progress);
        }
        if ( (iter+1) % 100 == 0 ) {
            for ( int i = 0; i < mod.n; ++i ) {
                if ( theta_accepts[i] < 20 ) {
                    mod.theta_sds[i] -= ((20.0 - theta_accepts[i]) * 0.01);
                }
                if ( theta_accepts[i] > 25 ) {
                    mod.theta_sds[i] += ((theta_accepts[i] - 25.0) * 0.01);
                }
                if ( mod.theta_sds[i] < 0 ) {
                    mod.theta_sds[i] = 0.01;
                }
                theta_accepts[i] = 0;
            }
            for ( int i = 0; i < mod.m; ++i ) {
                if ( alpha_accepts[i] < 20 ) {
                    mod.alpha_sds[i] -= ((20.0 - alpha_accepts[i]) * 0.01);
                }
                if ( alpha_accepts[i] > 25 ) {
                    mod.alpha_sds[i] += ((alpha_accepts[i] - 25.0) * 0.01);
                }
                if ( mod.alpha_sds[i] < 0 ) {
                    mod.alpha_sds[i] = 0.01;
                }
                alpha_accepts[i] = 0;
                if ( delta_accepts[i] < 20 ) {
                    mod.delta_sds[i] -= ((20.0 - delta_accepts[i]) * 0.01);
                }
                if ( delta_accepts[i] > 25 ) {
                    mod.delta_sds[i] += ((delta_accepts[i] - 25.0) * 0.01);
                }
                if ( mod.delta_sds[i] < 0 ) {
                    mod.delta_sds[i] = 0.01;
                }
                delta_accepts[i] = 0;
                int Kj = mod.K[i] - 1;
                if ( tau_accepts[i] < (20 * Kj) ) {
                    mod.tau_sds[i] -= (((20.0 * Kj) - tau_accepts[i])
                                       * (0.01 / Kj));
                }
                if ( tau_accepts[i] > (25 * Kj) ) {
                    mod.tau_sds[i] += ((tau_accepts[i] - (25.0 * Kj))
                                       * (0.01 / Kj));
                }
                if ( mod.tau_sds[i] < 0 ) {
                    mod.tau_sds[i] = 0.01;
                }
                tau_accepts[i] = 0;
            }
        }
        Rcpp::NumericVector old_theta = mod.theta(Rcpp::_, 0);
        mod.update_theta(0);
        for ( int i = 0; i < mod.n; ++i ) {
            if ( old_theta[i] != mod.theta(i, 0) ) {
                theta_accepts[i] += 1;
            }
        }
        Rcpp::NumericVector old_alpha = mod.alpha(Rcpp::_, 0);
        Rcpp::NumericVector old_delta = mod.delta(Rcpp::_, 0);
        ragged_array old_tau = mod.tau[0];
        mod.update_alpha(0);
        mod.update_delta(0);
        mod.update_tau(0);
        for ( int j = 0; j < mod.m; j++ ) {
            if ( old_alpha[j] != mod.alpha(j, 0) ) {
                alpha_accepts[j] += 1;
            }
            if ( old_delta[j] != mod.delta(j, 0) ) {
                delta_accepts[j] += 1;
            }
            for ( int k = 1; k < mod.K[j]; k++ ) {
                if ( old_tau[j][k] != mod.tau[0][j][k] ) {
                    tau_accepts[j] += 1;
                }
            }
        }
    }
    // Close out the progress display
    Rprintf("\rTuning proposals:    100.000 %%\n");
    // And return the results
    return Rcpp::List::create(mod.theta_sds,
                              mod.alpha_sds,
                              mod.delta_sds,
                              mod.tau_sds);
}
