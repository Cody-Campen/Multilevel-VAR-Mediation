#' simulate_data function
#' 
#' This function simulates an example dataset for mediation through person-specific VAR parameters.
#' 
#' @param n_people The number of people in the dataset. 
#' @param n_times The number of time points per person.
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
simulate_data = function(n_people = 100,
                         n_times = 30,
                         percent_missing = 0,
                         n_treatments = 2,
                         n_mediators = 2,
                         n_outcomes = 1,
                         treatment_effect_matrix = matrix(c(0, 0, 0,  # X -> M_1 intercept
                                                            0, 0, 0,  # X -> M_2 intercept
                                                            0, 0, 0,  # X -> M_1 autoregression
                                                            0, 0, 0,  # X -> M_1 to M_2 crossregression  
                                                            0, 0, 0,  # X -> M_2 to M_1 crossregression
                                                            0, 0, 0), # X -> M_2 autoregression 
                                                          nrow = 6, ncol = 3, byrow = T),
                         mediator_effect_matrix = matrix(c(0,   # M_1 intercept -> Y
                                                           0,   # M_2 intercept -> Y
                                                           0,   # M_1 autoregression -> Y
                                                           0,   # M_1 -> M_2 crossregression -> Y 
                                                           0,   # M_2 -> M_1 crossregression -> Y
                                                           0),  # M_2 autoregression -> Y
                                                       nrow = 6, ncol = 1, byrow = T),
                         direct_effect = matrix(c(0,  # Y intercept
                                                  0,  # X1 -> Y
                                                  0), # X2 -> Y
                                                nrow = 3, ncol = 1, byrow = T),
                         parameter_matrix_covariance = diag(6),
                         Y_covariance = diag(1),
                         m_standard_deviations = c(1, 1),
                         m_correlations = .2){
  # ---- Errors and warnings ----
  
  # Checking on the between_coefficient_matrix
  if(!is.matrix(treatment_effect_matrix)){
    stop("The argument for 'treatment_effect_matrix' is not a matrix")
  }
  
  if(!is.matrix(mediator_effect_matrix)){
    stop("The argument for 'mediator_effect_matrix' is not a matrix")
  }
  
  # ---- Generating the dataset ----
  
  id = rep(1:n_people, each = n_times)
  time = 1:n_times
  n_parameters = n_mediators + n_mediators^2
  
  X = mvrnorm(n = n_people, mu = rep(0, times = n_treatments), Sigma = diag(n_treatments))
  X = cbind(rep(1, n_people), X) # including the 1s to model the intercepts
  parameter_matrix_error = t(mvrnorm(n = n_people, mu = rep(0, times = n_parameters), Sigma = parameter_matrix_covariance))
    
  # First, define the between-level model, the effect of X on the person-specific parameters.
  parameter_matrix = treatment_effect_matrix %*% t(X) + parameter_matrix_error
  
  # Second, use these parameters to get the value of the mediator by...
  
  # ...extracting the intercepts and coefficients in the parameter matrix into more interpretable units.
  m_intercepts = parameter_matrix[1:n_mediators,]
  
  m_transition_matrix = array(parameter_matrix[(n_mediators+1):(n_mediators^2 + n_mediators),], 
                              dim = c(n_mediators, n_mediators, n_people))
  
  # and get the final Pearson correlations in the transition matrix from the fisher z transforms
  for(this_mediator in 1:n_mediators){
    m_transition_matrix[this_mediator, this_mediator, ] = tanh(m_transition_matrix[this_mediator, this_mediator, ])
  }
  
  # and of course, check for stability.
  for(this_person in 1:n_people) {
    eigen_values = eigen(m_transition_matrix[, , this_person])$values
    max_eigen_value = max(Mod(eigen_values))
    
    # when a participant has an eigen value greater than .99, this resamples their transition matrix until its below .99
    if(max_eigen_value >= .99){
      counter = 0
      while(max_eigen_value >= .99){
        # Let the user know when a transition matrix gets redrawn + how many times its getting redrawn
        counter = 1 + counter
        cat("\rRedrawing transition matrix for participant ",this_person, ". ", "Redraw number: ", counter, sep = "")
        
        # Resample the parametermatrix
        parameter_matrix_error[ , this_person] = t(mvrnorm(n = 1, mu = rep(0, times = n_parameters), Sigma = parameter_matrix_covariance))
        parameter_matrix[ , this_person] = treatment_effect_matrix %*% X[this_person, ] + parameter_matrix_error[ , this_person]
        m_transition_matrix[ , , this_person] = array(parameter_matrix[(n_mediators+1):(n_mediators^2 + n_mediators), this_person], 
                                                  dim = c(n_mediators, n_mediators))
        
        # Transform the value from normal to Pearson correlations
        for(this_mediator in 1:n_mediators){
          m_transition_matrix[this_mediator, this_mediator, this_person] = tanh(m_transition_matrix[this_mediator, this_mediator, this_person])
        }
        
        eigen_values = eigen(m_transition_matrix[, , this_person])$values
        max_eigen_value = max(Mod(eigen_values))
      }
      cat("\n")
    }
  }
  
  m_noise_covariance = array(NA, dim = c(n_mediators, n_mediators))
  m_process_noise_matrix = array(NA, dim = c(n_times, n_mediators))
  
################################################################################
#### This next chunk of code can ONLY accommodate 2 mediators at the moment ####
########### A change must be made to accommodate a different amount. ###########
################################################################################
  off_diagonal_index = !diag(nrow(m_noise_covariance))  
  diag(m_noise_covariance[ , ]) = m_standard_deviations^2
  m_noise_covariance[,][off_diagonal_index] =  prod(m_standard_deviations) * m_correlations
  m_process_noise_matrix[,] = mvrnorm(n = n_times, mu = rep(0, times = n_mediators), Sigma = m_noise_covariance)

  
  # Third, create a storage objects for the mediator time series...
  M = array(NA, dim = c(n_times, n_mediators, n_people)) # (timepoints x mediators x subjects)
  
  # ...and set the initial values for each subject.
  M[1, , ] = m_intercepts + m_process_noise_matrix[1,]
  
  # fourth, loop over timepoints to fill in participant's values for each their mediating time series
  for(this_person in 1:n_people){
    for(this_time in 2:n_times){
        M[this_time, ,this_person] = m_intercepts[, this_person] + 
                          m_transition_matrix[, , this_person] %*% (M[this_time - 1, ,this_person] - m_intercepts[, this_person]) + 
                          m_process_noise_matrix[this_time, ]
    }
  }
  
  # Finally, let's convert these transition matrix values to the outcome.
  Y_error = mvrnorm(n = n_people, mu = rep(0, times = n_outcomes), Sigma = Y_covariance)
  # add intercepts to the parameter matrix
  Y = t(parameter_matrix) %*% mediator_effect_matrix + X %*% direct_effect + Y_error

    # ---- Organizing and returning results ----
  dataset = data.frame(id = as.factor(id),
                       time = time)
  
  # Repeat X across each timepoint for long structure
  X_df = apply(X[,-1, drop = F], 2, function(x) rep(x, each = n_times))
  colnames(X_df) = paste("X", 1:n_treatments, sep = "")
  dataset = cbind(dataset, X_df)
  

  # Permute M to harmonize with the dataset structure
  M_df = matrix(aperm(M, c(1, 3, 2)), nrow = n_times*n_people, ncol = 2)
  
  # Create the missingness indices for M
  is_missing = as.logical(rbinom(n_people*n_times, 1, percent_missing))
  M_df[which(is_missing),] = NA
  
  colnames(M_df) = paste("M", 1:n_mediators, sep ="")
  dataset = cbind(dataset, M_df)
  
  # Repeat Y across timepoints for long structure
  Y_df = apply(Y, 2, function(y) rep(y, each = n_times))
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

