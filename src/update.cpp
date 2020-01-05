#include "model.h"

void Model::update_theta(int t) {
    for ( int i = 0; i < n; ++i ) {
        double proposal       = R::rnorm(theta(i, t), theta_sds[i]);
        double proposal_prior = theta_prior(proposal);
        double current_prior  = theta_prior(theta(i, t));
        double proposal_ll    = theta_ll(proposal, t, i);
        double current_ll     = theta_ll(theta(i, t), t, i);
        double accept_p       = proposal_ll - current_ll
                                + proposal_prior - current_prior;
        accept_p *= temps[t];
        if ( accept_p > 0.0 || log(R::runif(0.0, 1.0)) < accept_p) {
            theta(i, t) = proposal;
        }
    }
}

void Model::update_alpha(int t) {
    for ( int j = 0; j < m; ++j ) {
        double proposal       = R::rnorm(alpha(j, t), alpha_sds[j]);
        double proposal_prior = alpha_prior(proposal);
        if ( proposal_prior == R_NegInf ) {
            alpha_auto_rejects[j] += 1;
            return;
        }
        double current_prior = alpha_prior(alpha(j, t));
        double proposal_ll   = alpha_ll(proposal, t, j);
        double current_ll    = alpha_ll(alpha(j, t), t, j);
        double accept_p      = proposal_ll - current_ll
                               + proposal_prior - current_prior;
        accept_p *= temps[t];
        if ( accept_p > 0.0 || log(R::runif(0.0, 1.0)) < accept_p) {
            alpha(j, t) = proposal;
        }
    }
}

void Model::update_delta(int t) {
    for ( int j = 0; j < m; ++j ) {
        double proposal       = R::rnorm(delta(j, t), delta_sds[j]);
        double proposal_prior = delta_prior(proposal);
        if ( proposal_prior == R_NegInf ) {
            return;
        }
        double current_prior = delta_prior(delta(j, t));
        double proposal_ll   = delta_ll(proposal, t, j);
        double current_ll    = delta_ll(delta(j, t), t, j);
        double accept_p      = proposal_ll - current_ll
                               + proposal_prior - current_prior;
        accept_p *= temps[t];
        if ( accept_p > 0.0 || log(R::runif(0.0, 1.0)) < accept_p) {
            delta(j, t) = proposal;
        }
    }
}

void Model::update_tau(int t) {
    for ( int j = 0; j < m; ++j ) {
        for ( int k = 1; k < K[j]; ++k ) {
            Rcpp::NumericVector proposal = Rcpp::clone(tau[t][j]);
            proposal[k] = R::rnorm(tau[t][j][k], tau_sds[j]);
            double proposal_prior = tau_prior(proposal[k]);
            if ( proposal_prior == R_NegInf ) {
                continue;
            }
            double current_prior = tau_prior(tau[t][j][k]);
            double proposal_ll   = tau_ll(proposal, t, j);
            double current_ll    = tau_ll(tau[t][j], t, j);
            double accept_p      = proposal_ll - current_ll
                                   + proposal_prior - current_prior;
            accept_p *= temps[t];
            if ( accept_p > 0.0 || log(R::runif(0.0, 1.0)) < accept_p ) {
                tau[t][j] = proposal;
            }
        }
    }
}

