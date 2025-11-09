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
#' nsubs: the number of subjects
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
  # ---- (1a) level-2 likelihood functions ----
  for(this_sub in 1:nsubs){
      # First we start out for the likelihood for the person-specific VAR parameters as a function of their treatment, X.
      parameter_matrix.hat[1:n_parameters, this_sub] = X_fixed_effect[1:n_parameters, 1:n_treatments] %*% X[this_sub, 1:n_treatments] 
      parameter_matrix[1:n_parameters, this_sub] ~ dmnorm(parameter_matrix.hat[1:n_parameters, this_sub], parameter_matrix.precision)
      
      # Then we can get the likelihood for the outcome as a function of the person-specific VAR parameters
      Y.hat[1, this_sub] = M_fixed_effect[1:n_outcomes, 1:n_parameters] %*% parameter_matrix[1:n_parameters, this_sub] + 
                           direct_effect[1:n_outcomes, 1:n_treatments] %*% X[this_sub, 1:n_treatments]
      Y[this_sub, 1:n_outcomes] ~ dnorm(Y.hat[1:n_outcomes, this_sub], Y.precision) # for univariate outcomes
      
      # Y[this_sub, 1:n_outcomes] ~ dmnorm(Y.hat[1:n_outcomes, this_sub], Y.precision) # for multivariate outcomes
  }

  # and finally, calculate the indirect effects by cycling through each permutation of outcome, treatment, and parameter
  for(this_outcome in 1:n_outcomes){
    for(this_parameter in 1:n_parameters){
      for(this_treatment in 1:n_treatments){
        indirect_effect[this_outcome, this_parameter, this_treatment] = M_fixed_effect[this_outcome, this_parameter] * X_fixed_effect[this_parameter, this_treatment]
      }
    }
  }
  
  # ---- (1b) level-1 likelihood functions ----
  for(this_sub in 1:nsubs){
    # First, feed the intercepts from the parameter matrix into its own object
    M_intercept[this_sub, 1:n_mediators] = parameter_matrix[1:n_mediators, this_sub]

    # Then cycle the coefficients into their own object column-wise
    for(j in 1:n_mediators){ # column j
      for(i in 1:n_mediators){ # row i
        # The fancy indexing for parameter_matrix basically just counts up one
        # at a time through each parameter_matrix entry starting at n_mediators+1
        # (where the transition matrix coefficients should begin)
        M_transition_matrix[this_sub, i, j] = parameter_matrix[n_mediators+(j-1)*n_mediators+i, this_sub]
      }
    }

    # # Now we're going to begin on creating the process noise covariance matrix by first
    # # making the process noise, the correlations, then combining them into the covariance matrix.

    # Then feed the process noises into their own object 
    M_process_noise[this_sub, 1:n_mediators] = parameter_matrix[(n_mediators + (n_mediators*n_mediators) + 1):(2*n_mediators + (n_mediators*n_mediators)), this_sub]

    # Then feed the correlations between different mediators into its own object.
    # Note that this requires multiple mediators and should be removed for univariate mediation models
    for(j in 1:(n_mediators - 1)){ # column j
      for(i in (j+1):n_mediators){ # row i
        # The crazy indexing here does a similar thing as before, filling in the off-
        # diagonal elements column wise with the elements of parameter_matrix starting
        # at the value 2*n_mediators+n_mediators^2, which is where the correlations should
        # begin. In short, it takes the correlations from parameter matrix and reshapes it
        # into matrix of correlations.
        M_correlation_matrix[this_sub,i,j] = parameter_matrix[2*n_mediators + (n_mediators*n_mediators) + (j-1)*n_mediators - (j*j - j)/2 + (i-j), this_sub]
      }
    }

    for(this_var in 1:(n_mediators-1)){
      for(other_var in (this_var+1):n_mediators){
        this_covariance[this_sub, this_var, other_var] = M_process_noise[this_sub, other_var] * M_process_noise[this_sub, other_var] * M_correlation_matrix[this_sub, other_var, this_var]
        M_process_noise_covariance[this_sub, this_var, other_var] = this_covariance[this_sub, this_var, other_var]
        M_process_noise_covariance[this_sub, other_var, this_var] = this_covariance[this_sub, this_var, other_var]
      }
    }
  }

  
  # ---- (2a) level-2 likelihood priors ----
  
  for(this_parameter in 1:n_parameters){
    for(this_treatment in 1:n_treatments){
      X_fixed_effect[this_parameter, this_treatment] ~ dnorm(X_fixed_effect_mean[this_parameter, this_treatment], X_fixed_effect_precision[this_treatment, this_treatment])
    }
  }
  parameter_matrix.precision ~ dwish(parameter_rate_matrix, n_parameters+3)
  
  for(this_outcome in 1:n_outcomes){
    for(this_parameter in 1:n_parameters){
      M_fixed_effect[this_outcome, this_parameter] ~ dnorm(M_fixed_effect_mean[this_outcome, this_parameter], M_fixed_effect_precision[this_parameter, this_parameter])
    }
  }
  
  for(this_outcome in 1:n_outcomes){
    for(this_treatment in 1:n_treatments){
      direct_effect[this_outcome, this_treatment] ~ dnorm(direct_effect_mean[this_outcome, this_treatment], direct_effect_precision[this_treatment, this_treatment])
    }
  }
  
  Y.precision ~ dgamma(.1, .1) # for univariate outcomes
  # Y.precision ~ dwish(Y_rate_matrix, n_outcomes+3) # for multivariate outcomes
  
  # ---- (2b) level-1 likelihood priors
}