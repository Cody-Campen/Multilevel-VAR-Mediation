#' jags_model
#' 
#' This is the JAGS model to perform mediation analysis through person-specific VAR parameters.
#' Simply, it looks at the indirect effect of a time-invariant treatment on an outcome through 
#' person-specific VAR parameters, like a person-specific mediator intercept, autoregression, or cross-regression.
#' 
#' Author: Cody A. Campen
#' Last updated: 11/18/2025
#' For details see:
#' Insert our citation in apa format
#'   with a hanging indent
#'
#' ---- The data inputs ----
#' X_i: the treatment variable for person i
#' M_ij: a mediator time series for person i at time j
#' Y_i: the outcome variable for person i
#' n_people: the number of people/clusters in your dataset
#' n_seen: A vector of length n_people. The number of completed observations from each participant
#' n_miss: A vector of length n_people. The number of missing datapoints from each participant
#' 
#' 
#' n_treatments: the number of treatment variables used
#' n_mediators: the number of mediator time series used
#' n_parameters: the number of parameters estimated---which is equal to n_mediators + n_mediators^2.
#' n_outcomes: the number of outcome variables used
#' 
#' ---- Parameters ----
#' Level-2 parameters
#' X_to_intercepts: The fixed effects coefficient matrix of the treatment to the mediator intercepts
#' X_to_ARs: The fixed effects coefficient matrix of the treatment to the mediator autoregressive parameters in the transition matrix
#' X_to_CRs: The fixed effects coefficient matrix of the treatment to the mediator crossregressive parameters in the transition matrix
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
      M_intercept.hat[1:n_mediators, this_person] = X_to_intercept[1:n_mediators, 1:n_treatments] %*% X[this_person, 1:n_treatments]
      M_AR.hat[1:n_mediators, this_person] = X_to_AR[1:n_mediators, 1:n_treatments] %*% X[this_person, 1:n_treatments]
      M_CR.hat[1:n_mediators, this_person] = X_to_CR[1:n_mediators, 1:n_treatments] %*% X[this_person, 1:n_treatments]
 
      # Then sample from the vector of all parameters so that we can have covariances between intercepts, ARs and CRs
      parameter_matrix.hat[1:n_parameters, this_person] = c(M_intercept.hat[1:n_mediators, this_person],
                                                            M_AR.hat[1:n_mediators, this_person],
                                                            M_CR.hat[1:n_mediators, this_person])
      
      parameter_matrix[1:n_parameters, this_person] ~ dmnorm(parameter_matrix.hat[1:n_parameters, this_person], parameter_matrix.precision[1:n_parameters, 1:n_parameters])
      
      # Rename each parameter for easy tracking in the level-1 model.
      M_intercept[this_person, 1:n_mediators] = parameter_matrix[1:2, this_person]
      M_AR[this_person, 1:n_mediators] = parameter_matrix[3:4, this_person]
      M_CR[this_person, 1:n_mediators] = parameter_matrix[5:6, this_person]

      # Now we can get the likelihood for the outcome as a function of the person-specific VAR parameters
      Y.hat[1, this_person] = intercept_to_Y[1:n_outcomes, 1:n_mediators] %*% M_intercept[this_person, 1:n_mediators] + 
                              AR_to_Y[1:n_outcomes, 1:n_mediators] %*% M_AR[this_person, 1:n_mediators] +
                              CR_to_Y[1:n_outcomes, 1:n_mediators] %*% M_CR[this_person, 1:n_mediators] +
                              direct_effect[1:n_outcomes, 1:n_treatments] %*% X[this_person, 1:n_treatments]
      
      Y[this_person, 1:n_outcomes] ~ dnorm(Y.hat[1:n_outcomes, this_person], Y.precision)
  }
  
  # ---- (1.2) level-1 likelihood functions ----
  for(this_person in 1:n_people){
    # Begin by setting the initial values for each mediator time series
    M[this_person, 1, 1:n_mediators] ~ dmnorm(M_intercept[this_person, 1:n_mediators], M.precision[1:n_mediators, 1:n_mediators])

    # then define the likelihood for the rest of the time series for the non-missing values...
    for(this_time in times_seen[this_person, 1:n_seen[this_person]]){
      M.hat[this_person, this_time, 1:n_mediators] = M_intercept[this_person, 1:n_mediators] + M_transition_matrix[this_person, 1:n_mediators, 1:n_mediators] %*% (M[this_person, this_time-1, 1:n_mediators] - M_intercept[this_person, 1:n_mediators])
      M[this_person, this_time, 1:n_mediators] ~ dmnorm(M.hat[this_person, this_time, 1:n_mediators], M.precision[1:n_mediators, 1:n_mediators])
    }

    # ...and for the missing values.
    for(this_time in times_missed[this_person, 1:n_miss[this_person]]){
      M.hat[this_person, this_time, 1:n_mediators] = M_intercept[this_person, 1:n_mediators] + M_transition_matrix[this_person, 1:n_mediators, 1:n_mediators] %*% (M[this_person, this_time-1, 1:n_mediators] - M_intercept[this_person, 1:n_mediators])
      
      for(this_mediator in 1:n_mediators){
        M[this_person, this_time, this_mediator] ~ dnorm(M.hat[this_person, this_time, this_mediator], process_noise[this_mediator]^(-2))
      }
    }
  
    # ---- (1.2.1) Map the ARs and CRs from the level-2 parameter matrix onto the VAR transition matrix ----
    for(this_mediator in 1:n_mediators){
      M_transition_matrix[this_person, this_mediator, this_mediator] = M_AR[this_person, this_mediator]
    }
    M_transition_matrix[this_person, 2, 1] = M_CR[this_person, 1] # M1 -> M2
    M_transition_matrix[this_person, 1, 2] = M_CR[this_person, 2] # M2 -> M1
  }
  
  for(this_outcome in 1:n_outcomes){
    for(this_parameter in 1:n_mediators){
      for(this_treatment in 1:n_treatments){
        intercept_indirect_effect[this_outcome, this_parameter, this_treatment] = X_to_intercept[this_parameter, this_treatment] * intercept_to_Y[this_outcome, this_parameter]
        AR_indirect_effect[this_outcome, this_parameter, this_treatment] = X_to_AR[this_parameter, this_treatment] * AR_to_Y[this_outcome, this_parameter]
        CR_indirect_effect[this_outcome, this_parameter, this_treatment] = X_to_CR[this_parameter, this_treatment] * CR_to_Y[this_outcome, this_parameter]
      }
    }
  }
  
  # ---- (1.2.2) Define the person-invariance M precision matrix ----
  
  # Here, we're wanting to get all of the parts necessary (process noise, correlations and covariance matrix) to build the final process noise precision matrix
  
  for(this_mediator in 1:n_mediators){
    process_noise[this_mediator] = exp(log_process_noise[this_mediator])
    M_covariance_matrix[this_mediator, this_mediator] = process_noise[this_mediator] * process_noise[this_mediator]
  }
  
  for(this_mediator in 1:(n_mediators-1)){
    for(other_mediator in (this_mediator+1):n_mediators){
      correlation_matrix[this_mediator, other_mediator] = tanh(fisher_z[this_mediator, other_mediator])
      
      covariance_entry[this_mediator, other_mediator] = process_noise[this_mediator] * correlation_matrix[this_mediator, other_mediator] * process_noise[other_mediator]
      M_covariance_matrix[this_mediator, other_mediator] = covariance_entry[this_mediator, other_mediator]
      M_covariance_matrix[other_mediator, this_mediator] = covariance_entry[this_mediator, other_mediator]
    }
  }
  
  M.precision = inverse(M_covariance_matrix[1:n_mediators, 1:n_mediators])


  # ---- (2.1) level-2 likelihood priors ----
  
  # For the fixed effect of the treatment on the VAR parameters...
  for(this_parameter in 1:n_mediators){
    for(this_treatment in 1:n_treatments){
      X_to_intercept[this_parameter, this_treatment] ~ dnorm(0, .01)
      X_to_AR[this_parameter, this_treatment] ~ dnorm(0, 1)
      X_to_CR[this_parameter, this_treatment] ~ dnorm(0, 1)
    }
  }
  
  parameter_matrix.precision ~ dwish(parameter_matrix.rate[1:n_parameters, 1:n_parameters], 6)
  
  # ...and for the fixed effect of the treatment and parameters on the outcome.
  for(this_outcome in 1:n_outcomes){
    for(this_mediator in 1:n_mediators){
      intercept_to_Y[this_outcome, this_mediator] ~ dnorm(0, .01) 
      AR_to_Y[this_outcome, this_mediator] ~ dnorm(0, .01)
      CR_to_Y[this_outcome, this_mediator] ~ dnorm(0, .01)
    }
    for(this_treatment in 1:n_treatments){
      direct_effect[this_outcome, this_treatment] ~ dnorm(0, .01)
    }
  }

  Y.precision ~ dgamma(10, 10) # for univariate outcomes

  # ---- (2.2) level-1 likelihood priors ----

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