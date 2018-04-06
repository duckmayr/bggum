#' GGUM MCMC Sampler
#'
#' MCMC sampler for the generalized graded unfolding model (GGUM), utilizing
#' a Metropolis-Hastings algorithm
#'
#' \code{ggumMCMC} provides \code{R} implementation of an MCMC sampler for
#' the GGUM, based heavily on the algorithm given in de la Torre et al (2006);
#' though the package allows parameter estimation from \code{R},
#' the functions are actually written in \code{C++} to allow for reasonable
#' execution time.
#' 
#' Our sampler creates random initial values for the parameters of the model,
#' according to their prior distributions.
#' At each iteration, new parameter values are proposed
#' from a normal distribution with a mean of the current parameter value,
#' and the proposal is accepted probabilistically using a standard
#' Metropolis-Hastings acceptance ratio.
#' During burn-in, the standard deviation of the proposal densities
#' are tuned to ensure that the acceptance rate is neither too high nor too low
#' (we keep the acceptance rate between 0.2 and 0.25),
#' and the parameter draws are not stored.
#' Then the proposal densities are fixed and an additional number of draws
#' equal to \code{sample_iterations} are stored in a numeric matrix.
#'
#' @param data A numeric matrix giving the response by each
#'   respondent to each item
#' @param sample_iterations A numeric vector of length one;
#'   the number of iterations the sampler should store
#' @param burn_iterations A numeric vector of length one; the number of
#'   "burn-in" iterations to run, during which parameter draws are not
#'   stored. Currently, proposal densities are tuned during burn-in.
#'
#' @return A chain matrix; a numeric matrix with \code{sample_iterations} rows
#'   and one column for every parameter of the model, so that each element
#'   of the matrix gives the value of a parameter for a particular iteration
#'   of the MCMC algorithm.
#'
#' @seealso \code{\link{ggumProbability}}, \code{\link{ggumMC3}}
#'
#' @references Roberts, James S., John R. Donoghue, and James E. Laughlin.
#'   2000. ``A General Item Response Theory Model for Unfolding
#'   Unidimensional Polytomous Responses." \emph{Applied Psychological
#'   Measurement} 24(1): 3--32.
#' @references de la Torre, Jimmy, Stephen Stark, and Oleksandr S.
#'   Chernyshenko. 2006. ``Markov Chain Monte Carlo Estimation of Item
#'   Parameters for the Generalized Graded Unfolding Model." \emph{Applied
#'   Psychological Measurement} 30(3): 216--232.
#'   algorithm
#' @export
ggumMCMC <- function(data, sample_iterations, burn_iterations) {
    return(.ggumMCMC(data, sample_iterations, burn_iterations))
}
