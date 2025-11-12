#' create_inits function
#' 
#' This function creates the initial values from the priors and dataset information
#' 
#' @param X_fixed_effect_mean The prior means for the X fixed effects coefficient matrix
#' @param X_fixed_effect_covariance
#' @param parameter_rate_matrix
#' @param M_fixed_effect_mean
#' @param M_fixed_effect_covariance
#' @param direct_effect_mean
#' @param direct_effect_covariance
#' @param n_chains
#' @param n_treatments
#' @param n_parameters
#' @param n_mediators
#' @param n_outcomes
#' 
#' @return initial_values_list, a list containing n_chains lists with the following items:
#' \describe{
#' \item{X_fixed_effect} The initial value for each element of the X fixed effect coefficient matrix.
#' \item}parameter_matrix.precision} The initial value for each element of parameter matrix precision matrix.
#' \item{M_fixed_effect} The initial value for each element of the M fixed effect coefficient matrix.
#' \item{direct_effect} The initial value for each element of the direct effect matrix.
#' \item{Y.precision} The initial value for the precision of the outcome, Y.
#' \item{.RNG.name} The RNG passed to JAGS, set to Mersenne Twister.
#' \item{.RNG.seed} The RNG seed passed to JAGS. The seed is Different for each chain.
#' }
#' 
#' @export
create_inits = function(X_fixed_effect_mean,
                        X_fixed_effect_covariance,
                        parameter_rate_matrix,
                        M_fixed_effect_mean,
                        M_fixed_effect_covariance,
                        direct_effect_mean,
                        direct_effect_covariance,
                        n_chains, 
                        n_treatments, 
                        n_parameters,
                        n_mediators,
                        n_outcomes){
  initial_values_list = NULL
  
  for(this_chain in 1:n_chains){
    # creating storage objects for our initial values
    X_fixed_effect = matrix(NA, nrow = n_parameters, ncol = n_treatments+1)
    M_fixed_effect = matrix(NA, nrow = n_outcomes, ncol = n_parameters)
    direct_effect = matrix(NA, nrow = n_outcomes, ncol = n_treatments+1)
    
    log_process_noise = vector(length = n_mediators)
    fisher_z = matrix(NA, nrow = n_mediators, ncol = n_mediators)
    
    # sampling the initial values...
    
    # ...for the parameters in the level-2 model...
    for(this_parameter in 1:n_parameters){
      X_fixed_effect[this_parameter,] = mvrnorm(mu = X_fixed_effect_mean[this_parameter, ], Sigma = X_fixed_effect_covariance)
    }
    parameter_matrix.precision = rWishart(n = 1, df = n_parameters+3, Sigma = solve(parameter_rate_matrix)) |> drop()
    for(this_outcome in 1:n_outcomes){
      M_fixed_effect[this_outcome, ] = mvrnorm(mu = M_fixed_effect_mean[this_outcome, ], Sigma =M_fixed_effect_covariance)
      direct_effect[this_outcome, ] = mvrnorm(mu = direct_effect_mean[this_outcome, ], Sigma = direct_effect_covariance)
    }
    Y.precision = rgamma(n = 1, shape = .1, rate = .1)
    
    # ...and the process noise parameters in the level-1 model.
    log_process_noise = mvrnorm(mu = rep(0, times = n_mediators), Sigma = diag(n_mediators))
    
    for(this_mediator in 1:(n_mediators-1)){
      for(other_mediator in (this_mediator+1):n_mediators){
        fisher_z[this_mediator, other_mediator] = rnorm(1)
      }
    }
    
    # And finally, organize the samples in a list to pass onto the larger list of initial values
    list_to_add = list(X_fixed_effect = X_fixed_effect,
                       parameter_matrix.precision = parameter_matrix.precision,
                       M_fixed_effect = M_fixed_effect,
                       direct_effect = direct_effect,
                       Y.precision = Y.precision,
                       .RNG.name = "base::Mersenne-Twister",
                       .RNG.seed = 1000*this_chain)
    
    initial_values_list[[this_chain]] = list_to_add
  }
  return(initial_values_list)
}