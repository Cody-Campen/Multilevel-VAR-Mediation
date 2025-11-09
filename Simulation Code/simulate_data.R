#' simulate_data function
#' 
#' This function simulates an example dataset for mediation through person-specific VAR parameters.
#' 
#' @param nsubs The number of subjects.
#' @param ntimes The number of time points per subject.
#' @param nmediators The number of time-series mediators.
#' @param percent_missing The percentage of observations missing completely at random for each subject.
#' @param n_treatments The number of treatments for the model.
#' @param n_mediators The number of mediating time series for the model.
#' @param n_outcomes The number of outcomes for the model.
#' @param treatment_effects_matrix A (number of predictors)x(number of parameters) coefficient matrix for X, with the first row designated for the intercept.
#' @param mediator_effect_matrix A (number of outcomes)x(number of mediators) coefficient matrix for M parameters on Y, with no column designated for the intercept.
#' @param direct_effect A (number of predictors)x(number of outcomes) coefficient matrix for the direct effect of X on Y. 
#' 
#' @return
#' \describe{
#' \item{dataset} The simulated dataset from the given coefficient matrices
#' \item{natural_indirect_effect} The natural indirect effect through each path, calculated from the supplied matrices
#' \item{parameter_matrix} A (number of parameters)x(number of participants) matrix of each of the participants' parameter values
#' }
#' 
#' @export
simulate_data = function(nsubs = 100,
                         ntimes = 30,
                         percent_missing = 0,
                         n_treatments = 2,
                         n_mediators = 2,
                         n_outcomes = 2,
                         treatment_effect_matrix = matrix(c(0, 0, 0,  # X -> M_1 intercept
                                                            0, 0, 0,  # X -> M_2 intercept
                                                            0, 0, 0,  # X -> M_1 autoregression
                                                            0, 0, 0,  # X -> M_1 to M_2 crossregression  
                                                            0, 0, 0,  # X -> M_2 to M_1 crossregression
                                                            0, 0, 0,  # X -> M_2 autoregression 
                                                            0, 0, 0,  # X -> M_1 process noise sd
                                                            0, 0, 0,  # X -> M_2 process noise sd
                                                            0, 0, 0), # X -> M_1 & M_2 correlation
                                                          nrow = 9, ncol = 3, byrow = T),
                         mediator_effect_matrix = matrix(c(0,   # M_1 intercept -> Y
                                                           0,   # M_2 intercept -> Y
                                                           0,   # M_1 autoregression -> Y
                                                           0,   # M_1 -> M_2 crossregression -> Y 
                                                           0,   # M_2 -> M_1 crossregression -> Y
                                                           0,   # M_2 autoregression -> Y
                                                           0,   # M_1 log of process noise sd -> Y
                                                           0,   # M_2 log of process noise sd -> Y
                                                           0),  # fisher z transformed M_1 M_2 correlation -> Y
                                                       nrow = 9, ncol = 1, byrow = T),
                         direct_effect = matrix(c(0,  # Y intercept
                                                  0,  # X1 -> Y
                                                  0), # X2 -> Y
                                                nrow = 3, ncol = 1, byrow = T),
                         parameter_matrix_covariance = diag(9),
                         Y_covariance = diag(2)){
  # ---- Errors and warnings ----
  
  # Checking on the between_coefficient_matrix
  if(!is.matrix(treatment_effect_matrix)){
    stop("The argument for 'treatment_effect_matrix' is not a matrix")
  }
  
  if(!is.matrix(mediator_effect_matrix)){
    stop("The argument for 'mediator_effect_matrix' is not a matrix")
  }
  
  # ---- Generating the dataset ----
  
  id = rep(1:nsubs, each = ntimes)
  time = 1:ntimes
  n_parameters = dim(treatment_effect_matrix)[1]
  
  X = mvrnorm(n = nsubs, mu = rep(0, times = n_treatments), Sigma = diag(n_treatments))
  X = cbind(rep(1, nsubs), X) # including the 1s to model the intercepts
  parameter_matrix_error = t(mvrnorm(n = nsubs, mu = rep(0, times = n_parameters), Sigma = parameter_matrix_covariance))
    
  # First, define the between-level model, the effect of X on the person-specific parameters.
  parameter_matrix = treatment_effect_matrix %*% t(X) + parameter_matrix_error
  
  # Second, use these parameters to get the value of the mediator by...
  
  # ...extracting the coefficients in the parameter matrix into more interpretable units...
  m_intercepts = parameter_matrix[1:n_mediators,]
  
  m_transition_matrix = array(parameter_matrix[(n_mediators+1):(n_mediators^2 + n_mediators),], 
                              dim = c(n_mediators, n_mediators, nsubs))
  
  # and of course, check for stability
  for(this_sub in 1:nsubs) {
    eigen_values = eigen(m_transition_matrix[, , this_sub])$values
    max_eigen_value = max(Mod(eigen_values))
    
    if(max_eigen_value >= .99) {
      warning(paste("Participant", this_sub, "has unstable transition matrix. Max |eigenvalue| =", round(max_eigen_value, 3)))
    }
  }
  
  # ...and then transform the log standard deviations and fisher z transformed correlations into the 
  # process noise covariance matrix using exponential and hyperbolic tangent functions respectively.
  m_standard_deviations = parameter_matrix[(n_mediators^2 + n_mediators + 1):(n_mediators^2 + 2*n_mediators),] |> exp()
  m_correlations = parameter_matrix[(n_mediators^2 + 2*n_mediators+1):dim(parameter_matrix)[1],] |> tanh()
  
  m_noise_covariance = array(NA, dim = c(n_mediators, n_mediators, nsubs))
  m_process_noise_matrix = array(NA, dim = c(ntimes, n_mediators, nsubs))
  
################################################################################
###### This next for loop can ONLY accommodate 2 mediators at the moment  ###### 
########### A change must be made to accommodate a different amount. ########### 
################################################################################
  off_diagonal_index = !diag(nrow(m_noise_covariance))  
  for(this_sub in 1:nsubs){
    diag(m_noise_covariance[ , , this_sub]) = m_standard_deviations[,this_sub]^2
    m_noise_covariance[, , this_sub][off_diagonal_index] =  prod(m_standard_deviations[,this_sub]) * m_correlations[this_sub]
    m_process_noise_matrix[, ,this_sub] = mvrnorm(n = ntimes, mu = rep(0, times = n_mediators), Sigma = m_noise_covariance[,,this_sub])
  }

  
  # Third, create a storage objects for the mediator time series...
  M = array(NA, dim = c(ntimes, n_mediators, nsubs)) # (timepoints x mediators x subjects)
  
  # ...and set the initial values for each subject.
  M[1, , ] = m_intercepts + m_process_noise_matrix[1,,]
  
  # fourth, loop over timepoints to fill in participant's values for each their mediating time series
  for(this_sub in 1:nsubs){
    for(this_time in 2:ntimes){
        M[this_time, ,this_sub] = m_intercepts[, this_sub] + 
                          m_transition_matrix[, , this_sub] %*% (M[this_time - 1, ,this_sub] - m_intercepts[, this_sub]) + 
                          m_process_noise_matrix[this_time, ,this_sub]
    }
  }
  
  # Finally, let's convert these transition matrix values to the outcome.
  Y_error = mvrnorm(n = nsubs, mu = rep(0, times = n_outcomes), Sigma = Y_covariance)
  # add intercepts to the parameter matrix
  Y = t(parameter_matrix) %*% mediator_effect_matrix + X %*% direct_effect + Y_error

    # ---- Organizing and returning results ----
  dataset = data.frame(id = id,
                       time = time)
  
  # Repeat X across each timepoint for long structure
  X_df = apply(X[,-1], 2, function(x) rep(x, each = ntimes))
  colnames(X_df) = paste("X", 1:n_treatments, sep = "")
  dataset = cbind(dataset, X_df)
  

  # Permute M to harmonize with the dataset structure
  M_df = matrix(aperm(M, c(1, 3, 2)), nrow = ntimes*nsubs, ncol = 2)
  
  # Create the missingness indices for M
  is_missing = as.logical(rbinom(nsubs*ntimes, 1, percent_missing))
  M_df[which(is_missing),] = NA
  
  colnames(M_df) = paste("M", 1:n_mediators, sep ="")
  dataset = cbind(dataset, M_df)
  
  # Repeat Y across timepoints for long structure
  Y_df = apply(Y, 2, function(y) rep(y, each = ntimes))
  colnames(Y_df) = paste("Y", 1:n_outcomes, sep = "")
  dataset = cbind(dataset, Y_df)
  
  natural_indirect_effect = array(NA, dim = c(n_outcomes, n_parameters, n_treatments))
  for(this_outcome in 1:n_outcomes){
    for(this_parameter in 1:n_parameters){
      for(this_treatment in 1:n_treatments){
        natural_indirect_effect[this_outcome, this_parameter, this_treatment] = treatment_effect_matrix[this_parameter, this_treatment+1] * mediator_effect_matrix[this_parameter, this_outcome]
      }
    }
  }
  
  
  return(list(dataset = dataset,
              natural_indirect_effect = natural_indirect_effect,
              parameter_matrix = parameter_matrix))
}

