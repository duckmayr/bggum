#include "../inst/include/ggumR.h" 

using namespace Rcpp;

// I'll need to use this to shuffle the vector of chain indices
// (see http://gallery.rcpp.org/articles/stl-random-shuffle/)
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

//' GGUM MC3
//'
//' Metropolis Coupled Markov Chain Monte Carlo Sampling for the GGUM
//'
//' @param data A numeric matrix giving the individuals' responses
//' @param iters A vector of length one giving the number of iterations
//' @param N The number of chains
//' @param W The period by which to attempt chain swaps; e.g. if W = 100,
//'   a state swap will be proposed between two randomly selected chains
//'   every 100 iterations
//' @param Temps An optional provision of the temperatures for the chains;
//'   if not provided, each temperature T_t for t < N is given by
//'   T_{t+1} * (t + 1), and T_N = 1.
//'
//' @return A numeric matrix giving the parameter values at each iteration
//'   for the cold chain
//' @export
//[[Rcpp::export]]
NumericMatrix ggumMC3(IntegerMatrix data, int iters, int r_one, int r_two,
        int N, int W, Nullable<NumericVector> Temps = R_NilValue){
    // set up temperatures
    int howmanyswaps = 0, coldswaps = 0, coldattempts = 0;
    int one = 0, two = 1;
    IntegerVector chains = seq_len(N) - 1;
    NumericVector temps(N, 1.0);
    if ( Temps.isNull() ) {
        for ( int t = 1; t < N; ++t ) {
            temps[t] = 1.0 / (temps[t-1] * (t + 1));
        }
    }
    else {
        temps = as<NumericVector>(Temps);
    }
    // set up initial values
    int n = data.nrow();
    int m = data.ncol();
    IntegerVector K(m);
    for ( int j = 0; j < m; ++j ) {
        K[j] = unique(na_omit(data(_, j))).size();
    }
    NumericMatrix thetas(N, n);
    NumericMatrix alphas(N, m);
    NumericMatrix deltas(N, m);
    List taus(N);
    for ( int t = 0; t < N; ++t ) {
        thetas(t, _) = rnorm(n, 0.0, 1.0);
        alphas(t, _) = rep(1.0, m);
        deltas(t, _) = r4beta(m, 2, 2, -5, 5);
        List t_t(m);
        for ( int j = 0; j < m; ++j) {
            NumericVector tau_j(K[j]);
            for ( int k = 1; k < K[j]; ++k ) {
                tau_j[k] = r_4beta(2, 2, -6, 6);
            }
            t_t[j] = tau_j;
        }
        taus[t] = t_t;
    }
    // set up a matrix to store the cold chain values
    NumericMatrix result(iters, n+ 2*m + sum(K));
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
                th = updateTheta(th, data(i, _), a_t, d_t, t_t, temp);
                thetas(t, i) = th;
            }
            for ( int j = 0; j < m; ++j) { // for each item
                double a_j = a_t[j];
                double d_j = d_t[j];
                NumericVector t_j = as<NumericVector>(t_t[j]);
                IntegerVector answers = data(_, j);
                d_j = updateDelta(d_j, answers, th_t, a_j, t_j, temp);
                for ( int k = 1; k < K[j]; ++k ) {
                    t_j[k] = updateTau(k, answers, th_t, a_j, d_j, t_j, temp);
                }
                a_j = updateAlpha(a_j, answers, th_t, d_j, t_j, temp);
                deltas(t, j) = d_j;
                t_t[j] = t_j;
                alphas(t, j) = a_j;
            }
            taus[t] = t_t;
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
            P1 = sum(log(dnorm(th1, 0.0, 1.0)));
            P1 += sum(log(d4beta(a1, 1.5, 1.5, 0.25, 4)));
            P1 += sum(log(d4beta(d1, 2, 2, -5, 5)));
            List t1 = as<List>(taus[one]);
            P2 = sum(log(dnorm(th2, 0.0, 1.0)));
            P2 += sum(log(d4beta(a2, 1.5, 1.5, 0.25, 4)));
            P2 += sum(log(d4beta(d2, 2, 2, -5, 5)));
            List t2 = as<List>(taus[two]);
            IntegerVector answers(n);
            NumericVector t_1j(m), t_2j(m);
            for ( int j = 0; j < m; ++j ) {
                t_1j = as<NumericVector>(t1[j]);
                t_2j = as<NumericVector>(t2[j]);
                for ( int k = 1; k < t_1j.size(); ++k ) {
                    P1 += d_4beta(t_1j[k], 2, 2, -6, 6);
                    P2 += d_4beta(t_2j[k], 2, 2, -6, 6);
                }
                answers = data(_, j);
                L1 += sum(log(probCol(answers, th1, a1[j], d1[j], t_1j)));
                L2 += sum(log(probCol(answers, th2, a2[j], d2[j], t_2j)));
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
    Rcout << howmanyswaps << " successful swaps occurred.\n";
    Rcout << coldswaps << " were with the cold chain.\n";
    Rcout << "(Out of " << coldattempts << " attempted cold chain swaps.)\n";
    return result;
}
