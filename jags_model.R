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
      
      # So we get the intercepts
      M_intercept.hat[1:n_mediators, this_person] = X_to_intercepts[1:n_mediators, 1:n_treatments] %*% X[this_person, 1:n_treatments]
      for(this_mediator in 1:n_mediators){
        M_intercept[this_person, this_mediator] ~ dnorm(M_intercept.hat[this_mediator, this_person], 0.1) #diffuse for intercepts
      }
      
      # Then the ARs 
      M_AR.hat[1:n_mediators, this_person] = X_to_ARs[1:n_mediators, 1:n_treatments] %*% X[this_person, 1:n_treatments]
      for(this_AR in 1:n_mediators){
        M_AR[this_person, this_AR] ~ dnorm(M_AR.hat[this_ARm this_person], 2) T(-1, 1) # tighter for ARs
      }
      
      # and finally the CRs (hard coded for 2 mediators, for now)
      M_CR.hat[1:n_mediators, this_person] = X_to_CRs[1:n_mediators, 1:n_treatments] %*% X[this_person, 1:n_treatments]
      for(this_CR in 1:n_mediators){ # number of CRs = n_mediators, for this case where n_mediators = 2
        M_CR[this_person, this_CR] ~ dnorm(M_CR.hat[this_CR, this_person], 2) # tighter for CRs
      }
      
      # and combine everything into one person-specific vector so we can model M_fixed_effect as one vector in Y.hat
      parameter_matrix[1:n_parameters, this_person] = c(M_intercept[this_person, 1:n_mediators],
                                                        M_AR[this_person, 1:n_mediators],
                                                        M_CR[this_person, 1:n_mediators])
      
      # Then we can get the likelihood for the outcome as a function of the person-specific VAR parameters
      Y.hat[1, this_person] = M_fixed_effect[1:n_outcomes, 1:n_parameters] %*% parameter_matrix[1:n_parameters, this_person] + direct_effect[1:n_outcomes, 1:n_treatments] %*% X[this_person, 1:n_treatments]
      Y[this_person, 1:n_outcomes] ~ dnorm(Y.hat[1:n_outcomes, this_person], Y.precision)
  }
  
  # ---- (1.2) level-1 likelihood functions ----
  for(this_person in 1:n_people){
    # Begin by setting the initial values for each mediator time series
    M[this_person, 1, 1:n_mediators] ~ dmnorm(M_intercept[this_person, 1:n_mediators], M.precision[1:n_mediators, 1:n_mediators])

    # then define the likelihood for the rest of the time series for the observed values
    for(this_time in times_seen[this_person, 1:n_seen[this_person]]){
      M.hat[this_person, this_time, 1:n_mediators] = M_intercept[this_person, 1:n_mediators] + M_transition_matrix[this_person, 1:n_mediators, 1:n_mediators] %*% (M[this_person, this_time-1, 1:n_mediators] - M_intercept[this_person, 1:n_mediators])
      M[this_person, this_time, 1:n_mediators] ~ dmnorm(M.hat[this_person, this_time, 1:n_mediators], M.precision[1:n_mediators, 1:n_mediators])
    }

    # and for the missing values
    for(this_time in times_missed[this_person, 1:n_miss[this_person]]){
      M.hat[this_person, this_time, 1:n_mediators] = M_intercept[this_person, 1:n_mediators] + M_transition_matrix[this_person, 1:n_mediators, 1:n_mediators] %*% (M[this_person, this_time-1, 1:n_mediators] - M_intercept[this_person, 1:n_mediators])

      for(this_mediator in 1:n_mediators){
        M[this_person, this_time, this_mediator] ~ dnorm(M.hat[this_person, this_time, this_mediator], pow(process_noise[this_mediator], -2))
      }
    }
  
    # ---- (1.2.1) Take the entries from the parameter matrix (estimated above on line 34) into their own objects with interpretable names ----
    for(this_mediator in 1:n_mediators){
      M_transition_matrix[this_person, this_mediator, this_mediator] = M_AR[this_person, this_mediator]
    }
    M_transition_matrix[this_person, 1, 2] = M_CR[this_person, 1]
    M_transition_matrix[this_person, 2, 1] = M_CR[this_person, 2]
  }
  
  for(this_outcome in 1:n_outcomes){
    for(this_mediator in 1:n_mediators){
      for(this_treatment in 1:n_treatments){
        intercepts_indirect_effect[this_outcome, this_mediator, this_treatment] = M_fixed_effect[this_outcome, this_mediator] * X_to_intercepts[this_mediator, this_treatment]
        AR_indirect_effect[this_outcome, this_mediator, this_treatment] = M_fixed_effect[this_outcome, n_mediators + this_mediator] * X_to_ARs[this_mediator, this_treatment]
        CR_indirect_effect[this_outcome, this_mediator, this_treatment] = M_fixed_effect[this_outcome, 2 * n_mediators + this_mediator] * X_to_CRs[this_mediator, this_treatment]
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
  
  # For the parameters used to generate the mediator (X_fixed_effect and parameter_matrix.precision)
  for(this_mediator in 1:n_mediators){
    X_to_intercepts[this_mediator, 1:n_treatments] ~ dmnorm(c(0,0,0), X_to_intercepts.precision)
    
    X_to_ARs[this_mediator, 1:n_treatments] ~ dmnorm(c(0,0,0), X_to_ARs.precision)
                                                     
    X_to_CRs[this_mediator, 1:n_treatments] ~ dmnorm(c(0,0,0), X_to_CRs.precision)
  }

  parameter_matrix.precision[1:n_parameters, 1:n_parameters] ~ dwish(parameter_rate_matrix, n_parameters+3)
  
  # And for the parameters used to generate the outcome (M_fixed_effect, direct_effect, and Y.precision)
  for(this_outcome in 1:n_outcomes){
    M_fixed_effect[this_outcome, 1:n_parameters] ~ dmnorm(M_fixed_effect_mean[this_outcome, 1:n_parameters], M_fixed_effect_precision) 
  }
  
  for(this_outcome in 1:n_outcomes){ 
    direct_effect[this_outcome, 1:n_treatments] ~ dmnorm(direct_effect_mean[this_outcome, 1:n_treatments], direct_effect_precision)
  }
  
  Y.precision ~ dgamma(.1, .1) # for univariate outcomes

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