#' jags_model
#' 
#' This is the JAGS model to perform mediation analysis through person-specific VAR parameters.
#' Simply, it looks at the indirect effect of a time-invariant treatment on an outcome through 
#' person-specific VAR parameters, like a person-specific mediator intercept, autoregression, or cross-regression.
#' 
#' Author: Cody A. Campen
#' Last updated: 11/12/2025
#' For details see:
#' Insert our citation in apa format
#'   with a hanging indent
#'
#' ---- The data inputs ----
#' X_i: the treatment variable for person i
#' M_ij: a mediator time series for person i at time j
#' Y_i: the outcome variable for person i
#' n_people: the number of people/clusters in your dataset
#' n_times: the number of observations for each person/cluster
#' n_treatments: the number of treatment variables used
#' n_mediators: the number of mediator time series used
#' n_parameters: the number of parameters estimated---which is equal to n_mediators + n_mediators^2.
#' n_outcomes: the number of outcome variables used
#' X_fixed_effect_mean: A matrix with dimensions n_parameters x n_treatments. The prior mean matrix for the X_fixed_effect coefficient matrix
#' X_fixed_effect_cholesky: A matrix with dimensions n_treatment x n_treatments. The cholesky decomposition of the prior X_fixed_effects covariance matrix
#' M_fixed_effect_mean: A matrix with dimensions n_outcomes x n_parameters The prior mean matrix for the M_fixed_effect coefficient matrix
#' M_fixed_effect_cholesky: A matrix with dimensions n_parameters x n_parameters The cholesky decomposition of the prior M_fixed_effect covariance matrix
#' direct_effect_mean: A matrix with dimensions n_outcomes x n_treatments The prior mean matrix for the direct_effect coefficient matrix
#' direct_effect_cholesky: A matrix with dimensions n_treatment x n_treatment The cholesky decomposition of the prior direct_effect covariance matrix
#' parameter_rate_matrix: A matrix with dimensions n_parameters x n_parameters. The rate matrix for the Wishart prior of parameter_matrix.precision.
#' 
#' ---- Parameters ----
#' Level-2 parameters
#' X_fixed_effect: The fixed effects coefficient matrix of X on the parameter matrix mean
#' parameter_matrix: The parameter matrix for each participant. It contains the mediator intercepts and transition matrix parameters
#' M_fixed_effect: The fixed effects coefficient matrix of the parameter matrix on the outcome. 
#' direct_effect: The direct effect coefficient matrix of the treatment on the outcome.
#' indirect_effect: The product of each element of X_fixed_effect and M_fixed_effect, matched by the parameter it corresponds to. The resulting matrix is n_outcomes x n_parameters x n_treatments.
#' 
#' Level-1 parameters
#' M_intercept: The person-specific intercept for each mediating time series.
#' M_transition_matrix: The person-specific transition matrix for the mediators.
#' log_process_noise: The log of each mediators process noise. The log transform is so that it is easier to sample from the prior.
#' process_noise: The process noise of each mediator. Transformed from the log_process_noise through exponentiation. 
#' fisher_z: The fisher-Z transformed correlation matrix for the mediators.
#' correlation_matrix: The correlation matrix for the mediators. Transformed from fisher_z by taking the hyperbolic tangent.
#' M_covariance_matrix: The covariance matrix for the mediators.
#' M.precision: The precision matrix for the mediators. Transformed from M_covariance_matrix by taking the inverse.
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
    
    # ---- (1.2.1) Take the entries from the parameter matrix (estimated above on line 34) into their own objects with interpretable names ----
    
    # for the intercept entries of the parameter matrix (the first n_mediators)
    M_intercept[this_person, 1:n_mediators] = parameter_matrix[1:n_mediators, this_person]
    
    # and then coefficients entries (the remaining)
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
  # Here, we're wanting to get all of the parts necessary to build the process noise precision matrix (which gets used on line 53 and 58 in the M likelihood)
  
  # We begin by getting the process noise SDs and variances
  for(this_mediator in 1:n_mediators){
    # first we transform the log process noise to get the real process noise
    process_noise[this_mediator] = exp(log_process_noise[this_mediator])
    
    # and fill in the covariance diagonal while we're still in this loop
    M_covariance_matrix[this_mediator, this_mediator] = process_noise[this_mediator] * process_noise[this_mediator]
  }
  
  # then we move on to the off diagonal elements of the process noise covariance matrix
  for(this_mediator in 1:(n_mediators-1)){
    for(other_mediator in (this_mediator+1):n_mediators){
      # first we transform back the fisher_z prior to get the correlations between mediators
      correlation_matrix[this_mediator, other_mediator] = tanh(fisher_z[this_mediator, other_mediator])
      
      # and use those correlations and the process noise to get the covariance between each mediator
      covariance_entry[this_mediator, other_mediator] = process_noise[this_mediator] * correlation_matrix[this_mediator, other_mediator] * process_noise[other_mediator]
      M_covariance_matrix[this_mediator, other_mediator] = covariance_entry[this_mediator, other_mediator]
      M_covariance_matrix[other_mediator, this_mediator] = covariance_entry[this_mediator, other_mediator]
    }
  }
  
  # and finally, we invert our covariance matrix to get our desired precision matrix
  M.precision = inverse(M_covariance_matrix[1:n_mediators, 1:n_mediators])  


  # ---- (2.1) level-2 likelihood priors ----
  
  # For the parameters used to generate the mediator (X_fixed_effect and parameter_matrix.precision)
  for(this_parameter in 1:n_parameters){
    for(this_treatment in 1:n_treatments){
      X_fixed_effect_raw[this_parameter, this_treatment] ~ dnorm(0, 1)
    }
    X_fixed_effect[this_parameter, 1:n_treatments] = X_fixed_effect_mean[this_parameter, 1:n_treatments] + X_fixed_effect_cholesky[1:n_treatments, 1:n_treatments] %*% X_fixed_effect_raw[this_parameter, 1:n_treatments]
  }
  parameter_matrix.precision[1:n_parameters, 1:n_parameters] ~ dwish(parameter_rate_matrix, n_parameters+3)
  
  # And for the parameters used to generate the outcome (M_fixed_effect, direct_effect, and Y.precision)
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

  # ---- (2.2) level-1 likelihood priors ----
  # Because the intercept and transition matrix were handled at level-2, we only need to do the process noise priors
  
  # the log of the person-invariant process noise
  for(this_mediator in 1:n_mediators){
    log_process_noise[this_mediator] ~ dnorm(0, 1)
  }
  
  # the fisher z transformed correlations
  for(this_mediator in 1:(n_mediators-1)){
    for(other_mediator in (this_mediator+1):n_mediators){
      fisher_z[this_mediator, other_mediator] ~ dnorm(0, 1)
    }
  }
  
}