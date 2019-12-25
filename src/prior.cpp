#include "model.h"

double Model::theta_prior(double val) {
    return R::dnorm(val, 0.0, 1.0, 1);
}

double Model::alpha_prior(double val) {
    return d_4beta(val, alpha_parameters[0], alpha_parameters[1],
                   alpha_parameters[2], alpha_parameters[3], 1);
}

double Model::delta_prior(double val) {
    return d_4beta(val, delta_parameters[0], delta_parameters[1],
                   delta_parameters[2], delta_parameters[3], 1);
}

double Model::tau_prior(double val) {
    return d_4beta(val, tau_parameters[0], tau_parameters[1],
                   tau_parameters[2], tau_parameters[3], 1);
}

