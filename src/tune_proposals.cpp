#include "../inst/include/bggum.h"

using namespace Rcpp;

// [[Rcpp::export(.tune_proposals)]]
List tune_proposals(const IntegerMatrix& responseMatrix, NumericVector& thetas,
        NumericVector& alphas, NumericVector& deltas, List& taus,
        const IntegerVector& K, const int tune_iters, int n, int m,
        double th_prior_mean, double th_prior_sd, double a_shape1,
        double a_shape2, double a_a, double a_b, double d_shape1,
        double d_shape2, double d_a, double d_b, double t_shape1,
        double t_shape2, double t_a, double t_b) {
    // set up progress display
    Rcout.precision(1);
    double adv_sdtune_progress = 1000.0 / tune_iters;
    Rcout << "\rTuning proposals:    0%";
    int current_break = 1;
    // set up vectors for proposal SDs
    NumericVector theta_SD(n, 0.01);
    NumericVector alpha_SD(m, 0.01);
    NumericVector delta_SD(m, 0.01);
    NumericVector tau_SD(m, 0.01);
    // set up vectors to store number of acceptances per cycle
    IntegerVector theta_acceptances(n);
    IntegerVector alpha_acceptances(m);
    IntegerVector delta_acceptances(m);
    IntegerVector tau_acceptances(m);
    for ( int iter = 0; iter < tune_iters; ++iter ) {
        if ( (iter+1) % 10 == 0 ) {
            // Check for user interruption
            checkUserInterrupt();
            // Update progress display
            Rcout << "\rTuning proposals:    " << std::fixed << current_break * adv_sdtune_progress << "%";
            current_break += 1;
        }
        if ( (iter+1) % 100 == 0 ) {
            for ( int i = 0; i < n; ++i ) {
                if ( theta_acceptances[i] < 20 ) {
                    theta_SD[i] -= ((20.0 - theta_acceptances[i]) * 0.01);
                }
                if ( theta_acceptances[i] > 25 ) {
                    theta_SD[i] += ((theta_acceptances[i] - 25.0) * 0.01);
                }
                if ( theta_SD[i] < 0 ) {
                    theta_SD[i] = 0.01;
                }
                theta_acceptances[i] = 0;
            }
            for ( int i = 0; i < m; ++i ) {
                if ( alpha_acceptances[i] < 20 ) {
                    alpha_SD[i] -= ((20.0 - alpha_acceptances[i]) * 0.01);
                }
                if ( alpha_acceptances[i] > 25 ) {
                    alpha_SD[i] += ((alpha_acceptances[i] - 25.0) * 0.01);
                }
                if ( alpha_SD[i] < 0 ) {
                    alpha_SD[i] = 0.01;
                }
                alpha_acceptances[i] = 0;
            }
            for ( int i = 0; i < m; ++i ) {
                if ( delta_acceptances[i] < 20 ) {
                    delta_SD[i] -= ((20.0 - delta_acceptances[i]) * 0.01);
                }
                if ( delta_acceptances[i] > 25 ) {
                    delta_SD[i] += ((delta_acceptances[i] - 25.0) * 0.01);
                }
                if ( delta_SD[i] < 0 ) {
                    delta_SD[i] = 0.01;
                }
                delta_acceptances[i] = 0;
            }
            for ( int i = 0; i < m; ++i ) {
                int multiplier = K[i];
                if ( tau_acceptances[i] < (20 * multiplier) ) {
                    tau_SD[i] -= (((20.0 * multiplier) - tau_acceptances[i]) * (0.01 / multiplier));
                }
                if ( tau_acceptances[i] > (25 * multiplier) ) {
                    tau_SD[i] += ((tau_acceptances[i] - (25.0 * multiplier)) * (0.01 / multiplier));
                }
                if ( tau_SD[i] < 0 ) {
                    tau_SD[i] = 0.01;
                }
                tau_acceptances[i] = 0;
            }
        }
        for ( int i = 0; i < n; ++i ) {
            double theta = thetas[i];
            thetas[i] = update_theta_MCMC(responseMatrix(i, _), theta, alphas,
                    deltas, taus, theta_SD[i], th_prior_mean, th_prior_sd);
            if ( theta != thetas[i] ) {
                theta_acceptances[i] += 1;
            }
        }
        for ( int j = 0; j < m; j++ ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            alphas[j] = update_alpha_MCMC(responseMatrix(_, j), thetas, alpha,
                    delta, taus[j], alpha_SD[j], a_shape1, a_shape2, a_a, a_b);
            if ( alpha != alphas[j] ) {
                alpha_acceptances[j] += 1;
            }
        }
        for ( int j = 0; j < m; j++ ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            deltas[j] = update_delta_MCMC(responseMatrix(_, j), thetas, alpha,
                    delta, taus[j], delta_SD[j], d_shape1, d_shape2, d_a, d_b);
            if ( delta != deltas[j] ) {
                delta_acceptances[j] += 1;
            }
        }
        for ( int j = 0; j < m; j++ ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            NumericVector thisTau = taus[j];
            for ( int k = 1; k < K[j]; k++ ) {
                double thisTau_k = thisTau[k];
                thisTau[k] = update_tau_MCMC(k, responseMatrix(_, j), thetas,
                        alpha, delta, thisTau, tau_SD[j], t_shape1, t_shape2,
                        t_a, t_b);
                if ( thisTau_k != thisTau[k] ) {
                    tau_acceptances[j] += 1;
                }
            }
            taus[j] = thisTau;
        }
    }
    // Close out the progress display
    Rcout << "\rTuning proposals:    " << std::fixed << 100.0 << "%\n";
    // And return the results
    return List::create(theta_SD, alpha_SD, delta_SD, tau_SD);
}
