#' This is the code which runs the JAGS model and calculates the ESS and Rhat
#' Author: Cody A. Campen
#' 
#' Citation:
#' The citation of our article in apa format
#'    (with a hanging indent)
#'    

library(rjags)
library(mnormt)
library(MASS)
library(dplyr)
source("postcalc.R")
source("simulate_data.R")

the_seed = 1000

# Dataset generation variables
n_people = 100
n_times = 30
percent_missing = 0
n_treatments = 2
n_mediators = 2
n_parameters = 6
n_outcomes = 1
treatment_effect_matrix = matrix(c(0, 0, 2,  # X -> M_1 intercept
                                   0, 0, 0,  # X -> M_2 intercept
                                   0, 0, 0,  # X -> M_1 autoregression
                                   .2, 0, 0,  # X -> M_1 to M_2 crossregression  
                                   0, 0, 0,  # X -> M_2 to M_1 crossregression
                                   0, .3, 0), # X -> M_2 autoregression 
                                 nrow = n_parameters, ncol = n_treatments+1, byrow = T)
mediator_effect_matrix = matrix(c(3,   # M_1 intercept -> Y
                                  0,   # M_2 intercept -> Y
                                  0,   # M_1 autoregression -> Y
                                  0,   # M_1 -> M_2 crossregression -> Y 
                                  -2,   # M_2 -> M_1 crossregression -> Y
                                  0),  # M_2 autoregression -> Y
                                nrow = n_parameters, ncol = n_outcomes, byrow = T)
direct_effect = matrix(c(0,  # Y intercept
                         -5,  # X1 -> Y
                         5), # X2 -> Y
                       nrow = n_treatments + 1, ncol = n_outcomes, byrow = T)
parameter_matrix_covariance = diag(6)/1
Y_covariance = diag(n_outcomes)/1

# JAGS model variables
informative_priors = F  # T/F whether or not to use informative priors
n_chains = 2            # the number of chains used for the MCMC

# ---- Structure dataset ----
data_info = simulate_data(n_people = n_people, 
                          n_times = n_times, 
                          percent_missing = 0,
                          n_treatments = n_treatments, 
                          n_mediators = n_mediators, 
                          n_outcomes = n_outcomes,
                          treatment_effect_matrix = treatment_effect_matrix,
                          mediator_effect_matrix = mediator_effect_matrix,
                          direct_effect = direct_effect,
                          parameter_matrix_covariance = parameter_matrix_covariance,
                          Y_covariance = Y_covariance)

dataset = data_info$dataset
indirect_effect = data_info$natural_indirect_effect
parameter_matrix = data_info$parameter_matrix[1:n_parameters,]

# Structure the treatment
X = cbind(rep(1,n_people), as.matrix(dataset[dataset$time == 1,c("X1", "X2")], nrows=n_people))

# Structure the mediator
M_obs <- array(as.matrix(dataset[,c("M1", "M2")]), dim = c(n_times, n_people, n_mediators)) |> aperm(c(2,1,3))
# you want the dimensions to be structured: n_people, n_times, n_mediators

# Structure the outcome
Y = as.matrix(dataset[dataset$time == 1, c("Y1")], nrows = n_people)

# some priors

if(informative_priors){
  # Effect of X on the parameter matrix
  X_fixed_effect_mean = treatment_effect_matrix
  X_fixed_effect_covariance = diag(n_treatments+1)
  X_fixed_effect_cholesky = chol(X_fixed_effect_covariance)
  parameter_rate_matrix = solve(parameter_matrix_covariance)
  
  # Effect of the parameter matrix on the outcome
  M_fixed_effect_mean = t(mediator_effect_matrix)
  M_fixed_effect_covariance  = diag(n_parameters)
  M_fixed_effect_cholesky = chol(M_fixed_effect_covariance)
  Y_rate_matrix = solve(Y_covariance)
  
  # Effect of X on the outcome
  direct_effect_mean = t(direct_effect)
  direct_effect_covariance = diag(10, n_treatments+1)
  direct_effect_cholesky = chol(direct_effect_covariance)
} else{
  # Effect of X on the parameter matrix
  X_fixed_effect_mean = matrix(rep(0, times = n_parameters * (n_treatments+1)), nrow = n_parameters, ncol = n_treatments+1, byrow = T)
  X_fixed_effect_covariance = diag(10, n_treatments+1)
  X_fixed_effect_cholesky = chol(X_fixed_effect_covariance)
  parameter_rate_matrix = diag(10,n_parameters)
  
  # Effect of the parameter matrix on the outcome
  M_fixed_effect_mean = matrix(rep(0, times = n_parameters * n_outcomes), nrow = n_outcomes, ncol = n_parameters)
  M_fixed_effect_covariance  = diag(10, n_parameters)
  M_fixed_effect_cholesky = chol(M_fixed_effect_covariance)
  Y_rate_matrix = diag(10,n_outcomes)
  
  # Effect of X on the outcome
  direct_effect_mean = matrix(rep(0, times = n_treatments+1), nrow = n_outcomes, ncol = n_treatments+1)
  direct_effect_covariance = diag(10, n_treatments+1)
  direct_effect_cholesky = chol(direct_effect_covariance)
}

jags_data = list(X = X,
                 M = M_obs,
                 Y = Y,
                 n_people = n_people,
                 n_times = n_times,
                 n_treatments = n_treatments+1, # plus 1 because of the intercept term
                 n_parameters = n_parameters,
                 n_mediators = n_mediators,
                 n_outcomes = n_outcomes,
                 X_fixed_effect_mean = X_fixed_effect_mean,
                 X_fixed_effect_cholesky = X_fixed_effect_cholesky, # feeding in the cholesky decompositions instead of precision/covariance matrix for the non-centered parameterization
                 M_fixed_effect_mean = M_fixed_effect_mean,
                 M_fixed_effect_cholesky = M_fixed_effect_cholesky,
                 direct_effect_mean = direct_effect_mean,
                 direct_effect_cholesky = direct_effect_cholesky, 
                 parameter_rate_matrix = parameter_rate_matrix)

set.seed(the_seed)
create_inits = function(n_chains, 
                        n_treatments, 
                        n_parameters, 
                        n_outcomes ){
  initial_values_list = NULL
  
  for(this_chain in 1:n_chains){
    # creating storage objects for our initial values
    X_fixed_effect_raw = matrix(NA, nrow = n_parameters, ncol = n_treatments+1)
    M_fixed_effect_raw = matrix(NA, nrow = n_outcomes, ncol = n_parameters)
    direct_effect_raw = matrix(NA, nrow = n_outcomes, ncol = n_treatments+1)
    # sampling the initial values
    for(this_parameter in 1:n_parameters){
      X_fixed_effect_raw[this_parameter,] = mvrnorm(mu = rep(0, times = n_treatments+1), Sigma = diag(n_treatments+1))
    }
    parameter_matrix.precision = rWishart(n = 1, df = n_parameters+3, Sigma = solve(parameter_rate_matrix)) |> drop()
    for(this_outcome in 1:n_outcomes){
      M_fixed_effect_raw[this_outcome, ] = mvrnorm(mu = rep(0, times = n_parameters), Sigma = diag(n_parameters))
      direct_effect_raw[this_outcome, ] = mvrnorm(mu = rep(0, times = n_treatments+1), Sigma = diag(n_treatments+1))
    }
    Y.precision = rgamma(n = 1, shape = .1, rate = .1)
    
    # Organizing the samples in a list to pass onto the larger list of initial values
    list_to_add = list(X_fixed_effect_raw = X_fixed_effect_raw,
                       parameter_matrix.precision = parameter_matrix.precision,
                       M_fixed_effect_raw = M_fixed_effect_raw,
                       direct_effect_raw = direct_effect_raw,
                       Y.precision = Y.precision,
                       .RNG.name = "base::Mersenne-Twister",
                       .RNG.seed = 1000*this_chain)
    
    initial_values_list[[this_chain]] = list_to_add
  }
  
 return(initial_values_list)
}

inits_list = create_inits(n_chains, n_treatments, n_parameters, n_outcomes)
jagsModel = jags.model(file = "jags_model.R", 
                       data = jags_data, 
                       inits = inits_list, 
                       n.chains = n_chains, 
                       n.adapt = 1000) #n.adapt = 4000

update(jagsModel, n.iter = 1000)

# warm up/burn in period 2500

parameterlist = c("X_fixed_effect", "M_fixed_effect", "indirect_effect", "direct_effect")

codaSamples = coda.samples(jagsModel, variable.names = parameterlist, n.iter = 10000, thin = 1) 
save(codaSamples, file = "codaSamples.RData")

resulttable = zcalc(codaSamples)
resulttable

plot(codaSamples[[1]])
