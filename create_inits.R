#' create_inits function
#' 
#' This function creates the initial values from the priors and dataset information
#' 
#' @param n_chains
#' @param n_treatments
#' @param n_parameters
#' @param n_mediators
#' @param n_outcomes
#' 
#' @return initial_values_list, a list containing n_chains lists with the following items:
#' \describe{
#' \item{X_to_intercepts} 
#' \item{X_to_intercepts.precision}
#' \item{X_to_ARs}
#' \item{X_to_ARs.precision}
#' \item{X_to_CRs}
#' \item{X_to_CRs.precision}
#' \item{parameter_matrix.precision} The initial value for each element of parameter matrix precision matrix.
#' \item{M_fixed_effect} The initial value for each element of the M fixed effect coefficient matrix.
#' \item{direct_effect} The initial value for each element of the direct effect matrix.
#' \item{Y.precision} The initial value for the precision of the outcome, Y.
#' \item{.RNG.name} The RNG passed to JAGS, set to Mersenne Twister.
#' \item{.RNG.seed} The RNG seed passed to JAGS. The seed is Different for each chain.
#' }
#' 
#' @export
create_inits = function(n_chains, 
                        n_treatments, 
                        n_mediators,
                        n_outcomes,
                        the_seed){
  initial_values_list = NULL
  
  for(this_chain in 1:n_chains){
    # Sample the initial values for the parameters in the level-2 model...
    # ...for the regression model for M...
    X_to_intercept = matrix(runif(n_mediators * n_treatments, -2, 2),
                            nrow = n_mediators,
                            ncol = n_treatments)
    
    X_to_AR = matrix(runif(n_mediators * n_treatments, -1, 1),
                            nrow = n_mediators,
                            ncol = n_treatments)
    
    parameter_matrix.precision = diag(c(runif(n_mediators, 0, .4), runif(n_mediators, 8, 12)))
    
    # ...and the regression model for Y.
    intercept_to_Y = matrix(runif(n_outcomes * n_mediators, -2, 2),
                            nrow = n_outcomes,
                            ncol = n_mediators)
    
    AR_to_Y = matrix(runif(n_outcomes * n_mediators, -2, 2),
                           nrow = n_outcomes,
                           ncol = n_mediators)
    
    direct_effect = matrix(runif(n_outcomes * n_treatments, -2, 2),
                           nrow = n_outcomes,
                           ncol = n_treatments)

    Y.precision = runif(1, 0, 5)
    
    # Then, sample the initial value for the level-1 process noise variables.
    log_process_noise = runif(n_mediators, 0, 1)
  
    fisher_z = runif(1, -2, 2)
    
    # And finally, organize the samples in a list to pass onto the larger list of initial values
    list_to_add = list(X_to_intercept = X_to_intercept,
                       X_to_AR = X_to_AR,
                       intercept_to_Y = intercept_to_Y,
                       AR_to_Y = AR_to_Y,
                       direct_effect = direct_effect,
                       parameter_matrix.precision = parameter_matrix.precision,
                       Y.precision = Y.precision,
                       log_process_noise = log_process_noise,
                       fisher_z = fisher_z,
                       .RNG.name = "base::Mersenne-Twister",
                       .RNG.seed = 1000 * this_chain + the_seed) 

        initial_values_list[[this_chain]] = list_to_add
  }
  return(initial_values_list)
}