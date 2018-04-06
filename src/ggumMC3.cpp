#include "../inst/include/ggum.h" 

using namespace Rcpp;

//[[Rcpp::export(.ggumMC3)]]
NumericMatrix ggumMC3(IntegerMatrix data, int iters, int N, int W,
        NumericVector temps, NumericMatrix thetas, NumericMatrix alphas,
        NumericMatrix deltas, List taus, IntegerVector K, int n, int m){
    // set up some record-keeping variables
    int howmanyswaps = 0, coldswaps = 0, coldattempts = 0;
    int one = 0, two = 1;
    IntegerVector chains = seq_len(N) - 1;
    // set up a matrix to store the cold chain values
    NumericMatrix result(iters, n+ 2*m + sum(K));
    // set up progress display
    int progress = 0;
    int leftover = iters % 100;
    int adv_prog = (iters - leftover) / 100;
    Rcout << "\rRunning sampler: 0%";
    // now run the sampler
    for ( int iter = 0; iter < iters; ++iter ) { // for each iteration
        for ( int t = 0; t < N; ++t ) { // for each temperature
            NumericVector th_t = thetas(t, _);
            NumericVector a_t = alphas(t, _);
            NumericVector d_t = deltas(t, _);
            List t_t = as<List>(taus[t]);
            double temp = temps[t];
            for ( int i = 0; i < n; ++i ) { // for each respondent
                double th = th_t[i];
                th = update_theta_MC3(th, data(i, _), a_t, d_t, t_t, temp, 1.0);
                thetas(t, i) = th;
            }
            for ( int j = 0; j < m; ++j) { // for each item
                double a_j = a_t[j];
                double d_j = d_t[j];
                NumericVector t_j = as<NumericVector>(t_t[j]);
                IntegerVector answers = data(_, j);
                d_j = update_delta_MC3(d_j, answers, th_t, a_j, t_j, temp, 0.2);
                for ( int k = 1; k < K[j]; ++k ) {
                    t_j[k] = update_tau_MC3(k, answers, th_t, a_j, d_j, t_j,
                                            temp, 0.2);
                }
                a_j = update_alpha_MC3(a_j, answers, th_t, d_j, t_j, temp, 0.2);
                deltas(t, j) = d_j;
                t_t[j] = t_j;
                alphas(t, j) = a_j;
            }
            taus[t] = t_t;
        }
        if ( (iter+1) % adv_prog == 0 ) {
            if ( (iter+1) >= leftover ) {
                progress += 1;
                Rcout << "\rRunning sampler: " << progress << "%";
            }
        }
        if ( iter % W == 0 ) {
            checkUserInterrupt();
            if ( one == 0 || two == 0 ) {
                coldattempts += 1;
            }
            double r = 0.0, P1 = 0.0, P2 = 0.0, L1 = 0.0, L2 = 0.0;
            NumericVector th1 = thetas(one, _);
            NumericVector th2 = thetas(two, _);
            NumericVector a1 = alphas(one, _);
            NumericVector a2 = alphas(two, _);
            NumericVector d1 = deltas(one, _);
            NumericVector d2 = deltas(two, _);
            P1 = sum(dnorm(th1, 0.0, 1.0, true));
            P1 += sum(d4beta(a1, 1.5, 1.5, 0.25, 4, true));
            P1 += sum(d4beta(d1, 2, 2, -5, 5, true));
            List t1 = as<List>(taus[one]);
            P2 = sum(dnorm(th2, 0.0, 1.0, true));
            P2 += sum(d4beta(a2, 1.5, 1.5, 0.25, 4, true));
            P2 += sum(d4beta(d2, 2, 2, -5, 5, true));
            List t2 = as<List>(taus[two]);
            IntegerVector answers(n);
            NumericVector t_1j(m), t_2j(m);
            for ( int j = 0; j < m; ++j ) {
                t_1j = as<NumericVector>(t1[j]);
                t_2j = as<NumericVector>(t2[j]);
                for ( int k = 1; k < t_1j.size(); ++k ) {
                    P1 += d_4beta(t_1j[k], 2, 2, -6, 6, 1);
                    P2 += d_4beta(t_2j[k], 2, 2, -6, 6, 1);
                }
                answers = data(_, j);
                L1 += sum(log_probCol(answers, th1, a1[j], d1[j], t_1j));
                L2 += sum(log_probCol(answers, th2, a2[j], d2[j], t_2j));
            }
            double T = temps[one] - temps[two];
            double Y = L2 + P2 - L1 - P1;
            r = T * Y;
            if ( log(R::runif(0, 1)) < r ) {
                howmanyswaps += 1;
                if ( one == 0 || two == 0 ) {
                    coldswaps += 1;
                }
                NumericVector tmpThetas = thetas(one, _);
                NumericVector tmpAlphas = alphas(one, _);
                NumericVector tmpDeltas = deltas(one, _);
                List tmpTaus = as<List>(taus[one]);
                thetas(one, _) = thetas(two, _);
                alphas(one, _) = alphas(two, _);
                deltas(one, _) = deltas(two, _);
                taus[one] = taus[two];
                thetas(two, _) = tmpThetas;
                alphas(two, _) = tmpAlphas;
                deltas(two, _) = tmpDeltas;
                taus[two] = tmpTaus;
            }
            one = (one + 1) % (N - 1);
            two = (two % (N - 1)) + 1;
        }
        for ( int i = 0; i < n; ++i ) {
            result(iter, i) = thetas(0, i);
        }
        int K_ind = 0;
        List coldTau = as<List>(taus[0]);
        for ( int j = 0; j < m; ++j ) {
            result(iter, n+j) = alphas(0, j);
            result(iter, n+m+j) = deltas(0, j);
            NumericVector thisTau = as<NumericVector>(coldTau[j]);
            for ( int k = 1; k < K[j]; ++k ) {
                result(iter, n+2*m+k+K_ind) = thisTau[k];
            }
            K_ind += K[j];
        }
    }
    // return the cold chain
    Rcout << "\n" << howmanyswaps << " successful swaps occurred.\n";
    Rcout << coldswaps << " were with the cold chain.\n";
    Rcout << "(Out of " << coldattempts << " attempted cold chain swaps.)\n";
    return result;
}
