#' jags_model
#' 
#' This is the JAGS model to perform mediation analysis through person-specific VAR parameters.
#' Simply, it looks at the indirect effect of a time-invariant treatment on an outcome through 
#' person-specific VAR parameters, like a person's strength of autoregression or coupling between variables.
#' 
#' Author: Cody A. Campen
#' Last updated: 11/9/2025
#' For details see:
#' Insert our citation in apa format
#'   with a hanging indent
#'
#' ---- The data inputs ----
#' X_i: the treatment variable for person i
#' M_ij: a mediator time series for person i at time j
#' Y_i: the outcome variable for person i
#' 
#' ---- Auxiliary inputs ----
#' n_people: the number of subjects
#' n_treatments: the number of treatment variables
#' n_mediators: the number of mediator variables
#' n_parameters: the number of person-specific VAR parameters estimated
#' n_outcomes: the number of outcome variables
#' 
#' ---- Priors ----
#' X_fixed_effect_mean: The prior vector for the fixed effect of X on M
#' X_fixed_effect_variance: The prior variance-covariance matrix for the fixed effect of X on M
#' 
model {
  # ---- (1.1) level-2 likelihood functions ----
  for(this_person in 1:n_people){
      # First we start out for the likelihood for the person-specific VAR parameters as a function of their treatment, X.
      parameter_matrix.hat[1:n_parameters, this_person] = X_fixed_effect[1:n_parameters, 1:n_treatments] %*% X[this_person, 1:n_treatments]
      parameter_matrix[1:n_parameters, this_person] ~ dmnorm(parameter_matrix.hat[1:n_parameters, this_person], parameter_matrix.precision[1:n_parameters, 1:n_parameters])
      
      # Then we can get the likelihood for the outcome as a function of the person-specific VAR parameters
      Y.hat[1, this_person] = M_fixed_effect[1:n_outcomes, 1:n_parameters] %*% parameter_matrix[1:n_parameters, this_person] + direct_effect[1:n_outcomes, 1:n_treatments] %*% X[this_person, 1:n_treatments]
      Y[this_person, 1:n_outcomes] ~ dnorm(Y.hat[1:n_outcomes, this_person], Y.precision) # for univariate outcomes
  }

  # and finally, calculate the indirect effects by cycling through each permutation of outcome, treatment, and parameter
  for(this_outcome in 1:n_outcomes){
    for(this_parameter in 1:n_parameters){
      for(this_treatment in 1:n_treatments){
        indirect_effect[this_outcome, this_parameter, this_treatment] = M_fixed_effect[this_outcome, this_parameter] * X_fixed_effect[this_parameter, this_treatment]
      }
    }
  }
  
  # ---- (1.2) level-1 likelihood functions ----
  for(this_person in 1:n_people){
    # Begin by setting the initial values for each mediator time series
    M[this_person, 1, 1:n_mediators] ~ dmnorm(M_intercept[this_person, 1:n_mediators], M.precision[1:n_mediators, 1:n_mediators])

    # then define the likelihood for the rest of the time series
    for(this_time in 2:n_times){
      M.hat[this_person, this_time, 1:n_mediators] = M_intercept[this_person, 1:n_mediators] + M_transition_matrix[this_person, 1:n_mediators, 1:n_mediators] %*% ( M[this_person, this_time-1, 1:n_mediators] - M_intercept[this_person, 1:n_mediators] )
      M[this_person, this_time, 1:n_mediators] ~ dmnorm(M.hat[this_person, this_time, 1:n_mediators], M.precision[1:n_mediators, 1:n_mediators])
    }
    
    # ---- (1.2.1) Place the entries from the parameter matrix (estimated above) into their own objects ----
    # for the intercept
    M_intercept[this_person, 1:n_mediators] = parameter_matrix[1:n_mediators, this_person]

    # and then coefficients
    for(j in 1:n_mediators){ # column j
      for(i in 1:n_mediators){ # row i
        # The fancy indexing for parameter_matrix basically just counts up one
        # at a time through each parameter_matrix entry starting at n_mediators+1
        # and enters them column-wise into M_transition_matrix
        M_transition_matrix[this_person, i, j] = parameter_matrix[n_mediators+(j-1)*n_mediators+i, this_person]
      }
    }
  }
  
  # ---- (1.2.2) Define the person-invariant process noise structure ----
  for(this_mediator in 1:n_mediators){
    # transform the log process noise to get the real process noise
    process_noise[this_mediator] = exp(log_process_noise[this_mediator])
    
    # and fill in the covariance diagonal while we're here
    M_covariance_matrix[this_mediator, this_mediator] = process_noise[this_mediator] * process_noise[this_mediator]
    
    # and transform the fisher z transform back to real correlations
    correlation_matrix[this_mediator] = tanh(fisher_z[this_mediator])
  }
  
  for(this_mediator in 1:(n_mediators-1)){
    for(other_mediator in (this_mediator+1):n_mediators){
      covariance_entry[this_mediator, other_mediator] = process_noise[this_mediator] * correlation_matrix[this_mediator] * process_noise[other_mediator]
      M_covariance_matrix[this_mediator, other_mediator] = covariance_entry[this_mediator, other_mediator]
      M_covariance_matrix[other_mediator, this_mediator] = covariance_entry[this_mediator, other_mediator]
    }
  }
  M.precision = inverse(M_covariance_matrix[1:n_mediators, 1:n_mediators])  


  # ---- (2.1) level-2 likelihood priors ----
  
  for(this_parameter in 1:n_parameters){
    for(this_treatment in 1:n_treatments){
      X_fixed_effect_raw[this_parameter, this_treatment] ~ dnorm(0, 1)
    }
    X_fixed_effect[this_parameter, 1:n_treatments] = X_fixed_effect_mean[this_parameter, 1:n_treatments] + X_fixed_effect_cholesky[1:n_treatments, 1:n_treatments] %*% X_fixed_effect_raw[this_parameter, 1:n_treatments]
  }
  parameter_matrix.precision[1:n_parameters, 1:n_parameters] ~ dwish(parameter_rate_matrix, n_parameters+3)
  
  for(this_outcome in 1:n_outcomes){
    for(this_parameter in 1:n_parameters){
      M_fixed_effect_raw[this_outcome, this_parameter] ~ dnorm(0, 1)
    }
    M_fixed_effect[this_outcome, 1:n_parameters] = M_fixed_effect_mean[this_outcome, 1:n_parameters] + M_fixed_effect_cholesky[1:n_parameters, 1:n_parameters] %*% M_fixed_effect_raw[this_outcome, 1:n_parameters]
  }
  
  for(this_outcome in 1:n_outcomes){ 
    for(this_treatment in 1:n_treatments){
      direct_effect_raw[this_outcome, this_treatment] ~ dnorm(0, 1)
    }
    direct_effect[this_outcome, 1:n_treatments] = direct_effect_mean[this_outcome, 1:n_treatments] + direct_effect_cholesky[1:n_treatments, 1:n_treatments] %*% direct_effect_raw[this_outcome, 1:n_treatments]
  }
  
  Y.precision ~ dgamma(.1, .1) # for univariate outcomes
  # Y.precision ~ dwish(Y_rate_matrix, n_outcomes+3) # for multivariate outcomes
  
  # ---- (2.2) level-1 likelihood priors ----
  
  for(this_mediator in 1:n_mediators){
    # the log of the person-invariant process noise
    log_process_noise[this_mediator] ~ dnorm(0, 1)
    
    # the fisher z transformed correlations
    fisher_z[this_mediator] ~ dnorm(0, 1)
  }
  

}