#include "../inst/include/bggum.h" 

using namespace Rcpp;

// [[Rcpp::export(.tune_temperatures)]]
NumericVector tune_temps(IntegerMatrix data, int n_temps, int temp_tune_iters,
        int n_draws, int n, int m, IntegerVector K, List SDs,
        double th_prior_mean, double th_prior_sd, double a_shape1,
        double a_shape2, double a_a, double a_b, double d_shape1,
        double d_shape2, double d_a, double d_b, double t_shape1,
        double t_shape2, double t_a, double t_b){
    NumericVector temps(n_temps, 1.0);
    NumericVector tmp(n_draws);
    // set up initial values
    int N = 2;
    NumericMatrix thetas(N, n);
    NumericMatrix alphas(N, m);
    NumericMatrix deltas(N, m);
    List taus(N);
    for ( int t = 0; t < N; ++t ) {
        thetas(t, _) = rnorm(n, th_prior_mean, th_prior_sd);
        alphas(t, _) = rep(1.0, m);
        deltas(t, _) = r4beta(m, d_shape1, d_shape2, d_a, d_b);
        List t_t(m);
        for ( int j = 0; j < m; ++j) {
            NumericVector tau_j(K[j]);
            for ( int k = 1; k < K[j]; ++k ) {
                tau_j[k] = r_4beta(t_shape1, t_shape2, t_a, t_b);
            }
            t_t[j] = tau_j;
        }
        taus[t] = t_t;
    }
    NumericVector theta_SDs = as<NumericVector>(SDs[0]);
    NumericVector alpha_SDs = as<NumericVector>(SDs[1]);
    NumericVector delta_SDs = as<NumericVector>(SDs[2]);
    NumericVector tau_SDs = as<NumericVector>(SDs[3]);
    double rho1 = -2.0;
    double z = 1.0; // what Atchade et al. call n
    temps[1] = 1.0 / (1.0 + exp(rho1));
    int one = 0, two = 1;
    // set up progress display
    Rcout.precision(1);
    double adv_prog = 100.0 / ((temp_tune_iters / 100.0) * (n_temps - 1));
    int current_break = 1;
    Rcout << "\rTuning temperatures: 0%";
    double alph = 0.0;
    for ( int T = 1; T < n_temps; ++T ) {
        for ( int iter = 0; iter < temp_tune_iters; ++iter ) {
            if ( (iter+1) % 100 == 0 ) {
                // Every 100 iterations we check for a user interrupt
                // and update the progress bar
                checkUserInterrupt();
                double current_prog = current_break * adv_prog;
                Rcout << "\rTuning temperatures: " << std::fixed << current_prog << "%";
                current_break += 1;
            }
            for ( int t = one; t < two + 1; ++t ) { // for each temperature
                int ind = t % 2;
                NumericVector th_t = thetas(ind, _);
                NumericVector a_t = alphas(ind, _);
                NumericVector d_t = deltas(ind, _);
                List t_t = as<List>(taus[ind]);
                double temp = temps[t];
                for ( int i = 0; i < n; ++i ) { // for each respondent
                    double th = th_t[i];
                    th = update_theta_MC3(th, data(i, _), a_t, d_t, t_t, temp,
                            theta_SDs[i], th_prior_mean, th_prior_sd);
                    thetas(ind, i) = th;
                }
                for ( int j = 0; j < m; ++j) { // for each item
                    double a_j = a_t[j];
                    double d_j = d_t[j];
                    NumericVector t_j = as<NumericVector>(t_t[j]);
                    IntegerVector answers = data(_, j);
                    d_j = update_delta_MC3(d_j, answers, th_t, a_j, t_j, temp,
                            delta_SDs[j], d_shape1, d_shape2, d_a, d_b);
                    for ( int k = 1; k < K[j]; ++k ) {
                        t_j[k] = update_tau_MC3(k, answers, th_t, a_j, d_j,
                                t_j, temp, tau_SDs[j], t_shape1, t_shape2,
                                t_a, t_b);
                    }
                    a_j = update_alpha_MC3(a_j, answers, th_t, d_j, t_j, temp,
                            alpha_SDs[j], a_shape1, a_shape2, a_a, a_b);
                    deltas(ind, j) = d_j;
                    t_t[j] = t_j;
                    alphas(ind, j) = a_j;
                }
                taus[ind] = t_t;
            }
            alph = temps[one] - temps[two];
            double P1 = 0.0, P2 = 0.0, L1 = 0.0, L2 = 0.0;
            NumericVector th1 = thetas(one % 2, _);
            NumericVector th2 = thetas(two % 2, _);
            NumericVector a1 = alphas(one % 2, _);
            NumericVector a2 = alphas(two % 2, _);
            NumericVector d1 = deltas(one % 2, _);
            NumericVector d2 = deltas(two % 2, _);
            P1 = sum(dnorm(th1, th_prior_mean, th_prior_sd, true));
            P1 += sum(d4beta(a1, a_shape1, a_shape2, a_a, a_b, true));
            P1 += sum(d4beta(d1, d_shape1, d_shape2, d_a, d_b, true));
            List t1 = as<List>(taus[one % 2]);
            P2 = sum(dnorm(th2, th_prior_mean, th_prior_sd, true));
            P2 += sum(d4beta(a2, a_shape1, a_shape2, a_a, a_b, true));
            P2 += sum(d4beta(d2, d_shape1, d_shape2, d_a, d_b, true));
            List t2 = as<List>(taus[two % 2]);
            IntegerVector answers(n);
            NumericVector t_1j(m), t_2j(m);
            for ( int j = 0; j < m; ++j ) {
                t_1j = as<NumericVector>(t1[j]);
                t_2j = as<NumericVector>(t2[j]);
                for ( int k = 1; k < t_1j.size(); ++k ) {
                    P1 += d_4beta(t_1j[k], t_shape1, t_shape2, t_a, t_b, 1);
                    P2 += d_4beta(t_2j[k], t_shape1, t_shape2, t_a, t_b, 1);
                }
                answers = data(_, j);
                L1 += sum(log_probCol(answers, th1, a1[j], d1[j], t_1j));
                L2 += sum(log_probCol(answers, th2, a2[j], d2[j], t_2j));
            }
            alph *= (L2 + P2 - L1 - P1);
            alph = exp(alph);
            if ( alph > 1.0 ) {
                alph = 1.0;
            }
            rho1 += (1.0 / z) * (alph - 0.2338071);
            temps[two] = temps[one] / (1.0 + exp(rho1));
            if ( iter >= temp_tune_iters - n_draws ) {
                tmp[temp_tune_iters - iter - 1] = temps[two];
            }
        }
        one += 1;
        two += 1;
        z += 1.0;
        rho1 = -2.0;
        temps[T] = mean(tmp);
        if ( two < n_temps ) {
            temps[two] = temps[one] / (1.0 + exp(rho1));
        }
    }
    Rcout << "\rTuning temperatures: 100.0%\n";
    return temps;
}

