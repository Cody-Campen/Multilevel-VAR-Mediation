#' This is the code which runs the JAGS model and calculates the ESS and Rhat
#' Author: Cody A. Campen
#' 
#' Citation:
#' The citation of our article in apa format
#'    (with a hanging indent)
#'    

library(rjags)
library(MASS)
library(ks)
source("postcalc.R")
source("simulate_data.R")
source("create_inits.R")

the_seed = 1001
informative_priors = T  # T/F whether or not to use informative priors
diffuseness = 5
output_name = paste0("answers_", informative_priors, "_", diffuseness)

# Dataset generation variables
n_people = 138
n_times = 56
percent_missing = .3
n_treatments = 2
n_mediators = 2
n_parameters = 6
n_outcomes = 1
treatment_effect_matrix = matrix(c(0, 0, 2,   # X -> M_1 intercept
                                   0, 0, 0,   # X -> M_2 intercept
                                   0, 0, 0,   # X -> M_1 autoregression
                                   0, 0, 0,   # X -> M_1 to M_2 crossregression
                                   0, 0, 0,   # X -> M_2 to M_1 crossregression
                                   0, .3, 0), # X -> M_2 autoregression
                                 nrow = 6, ncol = n_treatments+1, byrow = T)
mediator_effect_matrix = matrix(c(3,   # M_1 intercept -> Y
                                  0,   # M_2 intercept -> Y
                                  0,   # M_1 autoregression -> Y
                                  0,   # M_1 -> M_2 crossregression -> Y
                                  -2,  # M_2 -> M_1 crossregression -> Y
                                  0),  # M_2 autoregression -> Y
                                nrow = 6, ncol = n_outcomes, byrow = T)
direct_effect = matrix(c(0,   # Y intercept
                         -3,  # X1 -> Y
                         3),  # X2 -> Y
                       nrow = n_treatments + 1, ncol = n_outcomes, byrow = T)
parameter_matrix_covariance = diag(6)/10 
Y_covariance = diag(n_outcomes)/1

# JAGS model variables
n_chains = 2

# ---- Structure dataset ----
set.seed(the_seed)
data_info = simulate_data(n_people = n_people, 
                          n_times = n_times, 
                          percent_missing = percent_missing,
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
M <- array(as.matrix(dataset[,c("M1", "M2")]), dim = c(n_times, n_people, n_mediators)) |> aperm(c(2,1,3))
# you want the dimensions to be structured: n_people, n_times, n_mediators

# Structure the outcome
Y = as.matrix(dataset[dataset$time == 1, c("Y1")], nrows = n_people)

# Now, we create objects to handle the time indexing
# Begin by creating objects to store our values
n_seen = vector(length = n_people)
n_miss = vector(length = n_people)
times_seen = matrix(NA, nrow = n_people, ncol = n_times)
times_missed = matrix(NA, nrow = n_people, ncol = n_times)

# and loop over each person to fill in their values
is_missing = matrix(NA, nrow = n_people, ncol = n_times)
for(this_person in 1:n_people){
  for(this_time in 2:n_times){
    # if any of the mediator values are missing at this time, then it will enter 1. If all are observed, it will enter 0.
    is_missing[this_person, this_time] = ifelse(sum(is.na(M[this_person, this_time, ]))==0, 0, 1)
  }
}

for(this_person in 1:n_people){
  n_miss[this_person] = length(which(is_missing[this_person,] %in% 1))
  # for those with missing values..
  if(n_miss[this_person] != 0){
    # record the location of their missing time point
    times_missed[this_person, 1:n_miss[this_person]] = which(is_missing[this_person,] %in% 1)
  }
  n_seen[this_person] = length(which(is_missing[this_person, ] %in% 0))
  times_seen[this_person, 1:n_seen[this_person]] = which(is_missing[this_person, ] %in% 0)
}


# ---- Define the priors ----

if(informative_priors){
  # Effect of X on the parameter matrix
  X_fixed_effect_mean = treatment_effect_matrix
  X_fixed_effect_covariance = diag(diffuseness, n_treatments+1)
  parameter_rate_matrix = diag(diffuseness, n_parameters)
  
  # Effect of the parameter matrix on the outcome
  M_fixed_effect_mean = t(mediator_effect_matrix)
  M_fixed_effect_covariance  = diag(diffuseness, n_parameters)
  Y_rate_matrix = diag(diffuseness, n_parameters)
  
  # Effect of X on the outcome
  direct_effect_mean = t(direct_effect)
  direct_effect_covariance = diag(diffuseness, n_treatments+1)
} else{
  # Effect of X on the parameter matrix
  X_fixed_effect_mean = matrix(rep(0, times = n_parameters * (n_treatments+1)), nrow = n_parameters, ncol = n_treatments+1, byrow = T)
  X_fixed_effect_covariance = diag(diffuseness, n_treatments+1)
  parameter_rate_matrix = diag(diffuseness, n_parameters)
  
  # Effect of the parameter matrix on the outcome
  M_fixed_effect_mean = matrix(rep(0, times = n_parameters * n_outcomes), nrow = n_outcomes, ncol = n_parameters)
  M_fixed_effect_covariance  = diag(diffuseness, n_parameters)
  Y_rate_matrix = diag(diffuseness,n_outcomes)
  
  # Effect of X on the outcome
  direct_effect_mean = matrix(rep(0, times = n_treatments+1), nrow = n_outcomes, ncol = n_treatments+1)
  direct_effect_covariance = diag(diffuseness, n_treatments+1)
}

# ---- Set up the initial values ----
inits_list = create_inits(X_fixed_effect_mean = X_fixed_effect_mean,
                          X_fixed_effect_covariance = X_fixed_effect_covariance,
                          parameter_rate_matrix = parameter_rate_matrix,
                          M_fixed_effect_mean = M_fixed_effect_mean,
                          M_fixed_effect_covariance = M_fixed_effect_covariance,
                          direct_effect_mean = direct_effect_mean,
                          direct_effect_covariance = direct_effect_covariance,
                          n_chains = n_chains, 
                          n_treatments = n_treatments, 
                          n_parameters = n_parameters,
                          n_mediators = n_mediators,
                          n_outcomes = n_outcomes)

jags_data = list(X = X,
                 M = M,
                 Y = Y,
                 n_people = n_people,
                 times_seen = times_seen,
                 times_missed = times_missed,
                 n_seen = n_seen,
                 n_miss = n_miss,
                 n_treatments = n_treatments+1, # plus 1 because of the intercept term
                 n_parameters = n_parameters,
                 n_mediators = n_mediators,
                 n_outcomes = n_outcomes,
                 X_fixed_effect_mean = X_fixed_effect_mean,
                 X_fixed_effect_precision = solve(X_fixed_effect_covariance), 
                 M_fixed_effect_mean = M_fixed_effect_mean,
                 M_fixed_effect_precision = solve(M_fixed_effect_covariance),
                 direct_effect_mean = direct_effect_mean,
                 direct_effect_precision = solve(direct_effect_covariance), 
                 parameter_rate_matrix = parameter_rate_matrix)


jagsModel = jags.model(file = "jags_model.R", 
                       data = jags_data, 
                       inits = inits_list, 
                       n.chains = n_chains, 
                       n.adapt = 4000)

update(jagsModel, n.iter = 15000)

# warm up/burn in period 2500

parameterlist = c("X_fixed_effect", "M_fixed_effect", "indirect_effect", "direct_effect")

before_time = Sys.time()
codaSamples = coda.samples(jagsModel, variable.names = parameterlist, n.iter = 20000, thin = 1) 
run_time = Sys.time() - before_time

resulttable = zcalc(codaSamples)

true_values = c(mediator_effect_matrix, vec(treatment_effect_matrix), direct_effect)
coda_answers = cbind(true_values = rep(true_values, length.out = nrow(resulttable)), resulttable)

answers = list(coda_answers = coda_answers,
               run_time = run_time)

assign(output_name, answers)

save(output_name, file = paste0("answers_","informative=",informative_priors, "_diffuseness=",diffuseness,".RData"))

# the dungeon

for(this_person in 1:n_people){
  for(this_time in times_missed[this_person, 1:n_miss[this_person]]){
    print(M[this_person, this_time-1, 1:n_mediators])
  }
}
