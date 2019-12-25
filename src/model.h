#ifndef BGGUM_MODEL_H
#define BGGUM_MODEL_H

#include <Rcpp.h>
#include <4beta.h>

// This is a helper function to extract a Rcpp::NumericVector from a Rcpp::List
inline Rcpp::NumericVector get_vec_from_list(Rcpp::List x, int i) {
    return Rcpp::as<Rcpp::NumericVector>(x[i]);
}

/* tau storage
 * We do not need to take advantage of the potential for data heterogeneity
 * that Rcpp::List was designed for; so, even though on the R side, tau will
 * be stored in a list, on the C++ side we will use ragged arrays.
 */
typedef std::vector<Rcpp::NumericVector> ragged_array;
typedef std::vector<std::vector<Rcpp::NumericVector> > nested_ragged_array;

class Model {

    public:

        // ATTRIBUTES
        //     - responses
        Rcpp::IntegerMatrix data;
        //     - parameters
        Rcpp::NumericMatrix theta;
        Rcpp::NumericMatrix alpha;
        Rcpp::NumericMatrix delta;
        nested_ragged_array tau;
        //     - size variables
        int n;
        int m;
        Rcpp::IntegerVector K;
        int T;
        //     - prior parameters
        Rcpp::NumericVector alpha_parameters;
        Rcpp::NumericVector delta_parameters;
        Rcpp::NumericVector tau_parameters;
        //     - proposal parameters
        Rcpp::NumericVector theta_sds;
        Rcpp::NumericVector alpha_sds;
        Rcpp::NumericVector delta_sds;
        Rcpp::NumericVector tau_sds;
        //     - NA index tracking
        std::vector<Rcpp::IntegerVector> item_non_na_indices;
        std::vector<Rcpp::IntegerVector> respondent_non_na_indices;
        //     - temperature schedule
        Rcpp::NumericVector temps;
        //     - state swap attempt tracking
        Rcpp::IntegerVector attempted_swaps;
        Rcpp::IntegerVector successful_swaps;
        
        // CONSTRUCTOR
        Model(Rcpp::IntegerMatrix data_,
              Rcpp::NumericMatrix theta_,
              Rcpp::NumericMatrix alpha_,
              Rcpp::NumericMatrix delta_,
              Rcpp::List tau_,
              Rcpp::IntegerVector K_,
              Rcpp::NumericVector theta_sds_,
              Rcpp::NumericVector alpha_sds_,
              Rcpp::NumericVector delta_sds_,
              Rcpp::NumericVector tau_sds_,
              Rcpp::NumericVector alpha_parameters_,
              Rcpp::NumericVector delta_parameters_,
              Rcpp::NumericVector tau_parameters_,
              Rcpp::NumericVector temps_) {

            // Set up the data and size variables
            data = data_;
            n = data.nrow();
            m = data.ncol();
            K = K_;
            T = temps_.size();

            // Set up the NA index tracking
            for ( int i = 0; i < n; ++i ) {
                Rcpp::IntegerVector these_indices;
                for ( int j = 0; j < m; ++j ) {
                    if ( !Rcpp::IntegerVector::is_na(data(i, j)) ) {
                        these_indices.push_back(j);
                    }
                }
                respondent_non_na_indices.push_back(these_indices);
            }
            for ( int j = 0; j < m; ++j ) {
                Rcpp::IntegerVector these_indices;
                for ( int i = 0; i < n; ++i ) {
                    if ( !Rcpp::IntegerVector::is_na(data(i, j)) ) {
                        these_indices.push_back(i);
                    }
                }
                item_non_na_indices.push_back(these_indices);
            }

            // Set up parameters and related variables
            theta = theta_;
            alpha = alpha_;
            delta = delta_;
            for ( int t = 0; t < tau_.size(); ++t ) {
                // Get the tau list for each parallel chain,
                // make it a vector of vectors,
                // then put them all together as a vector of vectors of vectors
                Rcpp::List l_tau_t = tau_[t];
                ragged_array tau_t;
                for ( int j = 0; j < m; ++j ) {
                    Rcpp::NumericVector tau_tj = get_vec_from_list(l_tau_t, j);
                    tau_t.push_back(tau_tj);
                }
                tau.push_back(tau_t);
            }
            alpha_parameters = alpha_parameters_;
            delta_parameters = delta_parameters_;
            tau_parameters = tau_parameters_;
            theta_sds = theta_sds_;
            alpha_sds = alpha_sds_;
            delta_sds = delta_sds_;
            tau_sds = tau_sds_;

            // Temperature schedule
            temps = temps_;

            // State swap attempt tracking
            attempted_swaps = Rcpp::rep(0, std::max(1, T-1));
            successful_swaps = Rcpp::rep(0, std::max(1, T-1));

        }

        // METHODS
        //     - (log) prior
        double theta_prior(double val);
        double alpha_prior(double val);
        double delta_prior(double val);
        double tau_prior(double val);
        //     - log likelihood
        double theta_ll(double val, int t, int i);
        double alpha_ll(double val, int t, int j);
        double delta_ll(double val, int t, int j);
        double tau_ll(Rcpp::NumericVector val, int t, int j);
        double item_ll(int t, int j);
        double log_likelihood(int t);
        //     - parameter updating
        void update_theta(int t);
        void update_alpha(int t);
        void update_delta(int t);
        void update_tau(int t);
        //     - attempt state swap
        void swap_states(int t1, int t2);

};

#endif

