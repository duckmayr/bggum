#include "model.h"

/*
 * The probability of a state swap between chains 1 and 2 is
 *  ( (prior1 * L1)^temp2 * (prior2 * L2)^temp1 )
 * -----------------------------------------------
 *  ( (prior1 * L1)^temp1 * (prior2 * L2)^temp2 )
 *
 * Then the log of this is
 *
 * temp2 * (log(prior1) + log(L1))
 * + temp1 * (log(prior2) + log(L2))
 * - temp1 * (log(prior1) + log(L1))
 * - temp2 * (log(prior2) + log(L2))
 *
 * which reduces to
 *
 * temp2 * (lp1 + ll1 - lp2 - ll2) + temp1 (lp2 + ll2 - lp1 - ll1)
 * = (temp1 - temp2) * (ll2 - ll1 + lp2 - lp1)
 */

void Model::swap_states_burn(int t1, int t2) {
    double chain1_prior = 0.0, chain2_prior = 0.0;
    for ( int i = 0; i < n; ++i ) {
        chain1_prior += theta_prior(theta(i, t1));
        chain2_prior += theta_prior(theta(i, t2));
    }
    for ( int j = 0; j < m; ++j ) {
        chain1_prior += alpha_prior(alpha(j, t1));
        chain1_prior += delta_prior(delta(j, t1));
        chain2_prior += alpha_prior(alpha(j, t2));
        chain2_prior += delta_prior(delta(j, t2));
        for ( int k = 1; k < K[j]; ++k ) {
            chain1_prior += tau_prior(tau[t1][j][k]);
            chain2_prior += tau_prior(tau[t2][j][k]);
        }
    }
    double chain1_log_likelihood = log_likelihood(t1);
    double chain2_log_likelihood = log_likelihood(t2);
    double accept_p = chain2_log_likelihood - chain1_log_likelihood;
    accept_p += (chain2_prior - chain1_prior);
    accept_p *= (temps[t1] - temps[t2]);
    if ( accept_p > 0.0 || log(R::runif(0.0, 1.0)) < accept_p ) {
        Rcpp::NumericVector tmp_theta = theta(Rcpp::_, t1);
        Rcpp::NumericVector tmp_alpha = alpha(Rcpp::_, t1);
        Rcpp::NumericVector tmp_delta = delta(Rcpp::_, t1);
        RaggedArray tmp_tau = tau[t1];
        theta(Rcpp::_, t1) = theta(Rcpp::_, t2);
        alpha(Rcpp::_, t1) = alpha(Rcpp::_, t2);
        delta(Rcpp::_, t1) = delta(Rcpp::_, t2);
        tau[t1] = tau[t2];
        theta(Rcpp::_, t2) = tmp_theta;
        alpha(Rcpp::_, t2) = tmp_alpha;
        delta(Rcpp::_, t2) = tmp_delta;
        tau[t2] = tmp_tau;
    }
}

void Model::swap_states(int t1, int t2) {
    attempted_swaps[t1] += 1;
    double chain1_prior = 0.0, chain2_prior = 0.0;
    for ( int i = 0; i < n; ++i ) {
        chain1_prior += theta_prior(theta(i, t1));
        chain2_prior += theta_prior(theta(i, t2));
    }
    for ( int j = 0; j < m; ++j ) {
        chain1_prior += alpha_prior(alpha(j, t1));
        chain1_prior += delta_prior(delta(j, t1));
        chain2_prior += alpha_prior(alpha(j, t2));
        chain2_prior += delta_prior(delta(j, t2));
        for ( int k = 1; k < K[j]; ++k ) {
            chain1_prior += tau_prior(tau[t1][j][k]);
            chain2_prior += tau_prior(tau[t2][j][k]);
        }
    }
    double chain1_log_likelihood = log_likelihood(t1);
    double chain2_log_likelihood = log_likelihood(t2);
    double accept_p = chain2_log_likelihood - chain1_log_likelihood;
    accept_p += (chain2_prior - chain1_prior);
    accept_p *= (temps[t1] - temps[t2]);
    if ( accept_p > 0.0 || log(R::runif(0.0, 1.0)) < accept_p ) {
        successful_swaps[t1] += 1;
        Rcpp::NumericVector tmp_theta = theta(Rcpp::_, t1);
        Rcpp::NumericVector tmp_alpha = alpha(Rcpp::_, t1);
        Rcpp::NumericVector tmp_delta = delta(Rcpp::_, t1);
        RaggedArray tmp_tau = clone(tau[t1]);
        theta(Rcpp::_, t1) = theta(Rcpp::_, t2);
        alpha(Rcpp::_, t1) = alpha(Rcpp::_, t2);
        delta(Rcpp::_, t1) = delta(Rcpp::_, t2);
        tau[t1] = tau[t2];
        theta(Rcpp::_, t2) = tmp_theta;
        alpha(Rcpp::_, t2) = tmp_alpha;
        delta(Rcpp::_, t2) = tmp_delta;
        tau[t2] = tmp_tau;
    }
}

