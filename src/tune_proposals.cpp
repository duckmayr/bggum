#include "../inst/include/ggum.h"

using namespace Rcpp;

//[[Rcpp::export]]
List tune_proposals(const IntegerMatrix& responseMatrix, NumericVector& thetas,
        NumericVector& alphas, NumericVector& deltas, List& taus,
        const IntegerVector& K, const int burn_iters, int n, int m) {
    // set up progress display
    Rcout.precision(1);
    double adv_sdtune_progress = 1000.0 / burn_iters;
    Rcout << "\rTuning proposals: 0%";
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
    for ( int iter = 0; iter < burn_iters; ++iter ) {
        if ( (iter+1) % 10 == 0 ) {
            // Check for user interruption
            checkUserInterrupt();
            // Update progress display
            Rcout << "\rTuning proposals: " << std::fixed << current_break * adv_sdtune_progress << "%";
            current_break += 1;
        }
        if ( (iter+1) % 100 == 0 ) {
            for ( int i = 0; i < n; ++i ) {
                if ( theta_acceptances[i] < 23 ) {
                    theta_SD[i] -= ((23.0 - theta_acceptances[i]) * 0.01);
                }
                if ( theta_acceptances[i] > 24 ) {
                    theta_SD[i] += ((theta_acceptances[i] - 25.0) * 0.01);
                }
                if ( theta_SD[i] < 0 ) {
                    theta_SD[i] = 0.01;
                }
                theta_acceptances[i] = 0;
            }
            for ( int i = 0; i < m; ++i ) {
                if ( alpha_acceptances[i] < 23 ) {
                    alpha_SD[i] -= ((23.0 - alpha_acceptances[i]) * 0.01);
                }
                if ( alpha_acceptances[i] > 24 ) {
                    alpha_SD[i] += ((alpha_acceptances[i] - 25.0) * 0.01);
                }
                if ( alpha_SD[i] < 0 ) {
                    alpha_SD[i] = 0.01;
                }
                alpha_acceptances[i] = 0;
            }
            for ( int i = 0; i < m; ++i ) {
                if ( delta_acceptances[i] < 23 ) {
                    delta_SD[i] -= ((23.0 - delta_acceptances[i]) * 0.01);
                }
                if ( delta_acceptances[i] > 24 ) {
                    delta_SD[i] += ((delta_acceptances[i] - 25.0) * 0.01);
                }
                if ( delta_SD[i] < 0 ) {
                    delta_SD[i] = 0.01;
                }
                delta_acceptances[i] = 0;
            }
            for ( int i = 0; i < m; ++i ) {
                int multiplier = K[i];
                if ( tau_acceptances[i] < (23 * multiplier) ) {
                    tau_SD[i] -= (((23.0 * multiplier) - tau_acceptances[i]) * (0.01 / multiplier));
                }
                if ( tau_acceptances[i] > (24 * multiplier) ) {
                    tau_SD[i] += ((tau_acceptances[i] - (24.0 * multiplier)) * (0.01 / multiplier));
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
                    deltas, taus, theta_SD[i]);
            if ( theta != thetas[i] ) {
                theta_acceptances[i] += 1;
            }
        }
        for ( int j = 0; j < m; j++ ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            alphas[j] = update_alpha_MCMC(responseMatrix(_, j), thetas, alpha,
                    delta, taus[j], alpha_SD[j]);
            if ( alpha != alphas[j] ) {
                alpha_acceptances[j] += 1;
            }
        }
        for ( int j = 0; j < m; j++ ) {
            double alpha = alphas[j];
            double delta = deltas[j];
            deltas[j] = update_delta_MCMC(responseMatrix(_, j), thetas, alpha,
                    delta, taus[j], delta_SD[j]);
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
                        alpha, delta, thisTau, tau_SD[j]);
                if ( thisTau_k != thisTau[k] ) {
                    tau_acceptances[j] += 1;
                }
            }
            taus[j] = thisTau;
        }
    }
    // Close out the progress display
    Rcout << "\rTuning proposals: " << std::fixed << 100.0 << "%\n";
    // And return the results
    return List::create(theta_SD, alpha_SD, delta_SD, tau_SD);
}
