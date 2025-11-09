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

the_seed = 1913
nsubs = 100
ntimes = 30
percent_missing = 0
n_treatments = 2
n_parameters = 9
n_mediators = 2
n_outcomes = 1
treatment_effect_matrix = matrix(c(0, 0, 0,  # X -> M_1 intercept
                                   0, 0, 0,  # X -> M_2 intercept
                                   0, 0, 5,  # X -> M_1 autoregression
                                   0, 0, 0,  # X -> M_1 to M_2 crossregression  
                                   0, 0, 0,  # X -> M_2 to M_1 crossregression
                                   0, 0, 0,  # X -> M_2 autoregression 
                                   0, 0, 0,  # X -> M_1 process noise sd
                                   0, 0, 0,  # X -> M_2 process noise sd
                                   0, 0, 0), # X -> M_1 & M_2 correlation
                                 nrow = n_parameters, ncol = n_treatments+1, byrow = T)
mediator_effect_matrix = matrix(c(0,   # M_1 intercept -> Y
                                  0,   # M_2 intercept -> Y
                                  3,   # M_1 autoregression -> Y
                                  0,   # M_1 -> M_2 crossregression -> Y 
                                  0,   # M_2 -> M_1 crossregression -> Y
                                  0,   # M_2 autoregression -> Y
                                  0,   # M_1 log of process noise sd -> Y
                                  0,   # M_2 log of process noise sd -> Y
                                  0),  # fisher z transformed M_1 M_2 correlation -> Y
                                nrow = n_parameters, ncol = n_outcomes, byrow = T)
direct_effect = matrix(c(0,  # Y intercept
                         5,  # X1 -> Y
                         0), # X2 -> Y
                       nrow = n_treatments + 1, ncol = n_outcomes, byrow = T)
parameter_matrix_covariance = diag(n_parameters)/100
Y_covariance = diag(n_outcomes)/100
# ---- Structure dataset ----

# For our dataset for...
data_info = simulate_data(nsubs = nsubs, ntimes = ntimes, percent_missing = 0,
                          n_treatments = n_treatments, n_outcomes = n_outcomes,
                          treatment_effect_matrix = treatment_effect_matrix,
                          mediator_effect_matrix = mediator_effect_matrix,
                          direct_effect = direct_effect,
                          parameter_matrix_covariance = parameter_matrix_covariance,
                          Y_covariance = Y_covariance)
dataset = data_info$dataset

indirect_effect = data_info$natural_indirect_effect
parameter_matrix = data_info$parameter_matrix[1:n_parameters,]

# Structure the treatment
X = cbind(rep(1,nsubs), as.matrix(dataset[dataset$time == 1,c("X1", "X2")], nrows=nsubs))

# Structure the mediator
M_obs <- array(as.matrix(dataset[,c("M1", "M2")]), dim = c(ntimes, nsubs, n_mediators)) |> aperm(c(2,1,3))
# you want the dimensions to be structured: nsubs, ntimes, n_mediators

# Structure the outcome
Y = as.matrix(dataset[dataset$time == 1, c("Y1")], nrows = nsubs)

# some priors

# ---- Effect of X on M mean and precision ----
# immensely informative priors
X_fixed_effect_mean = treatment_effect_matrix
X_fixed_effect_precision = diag(n_treatments+1)
parameter_rate_matrix = solve(parameter_matrix_covariance)

# # an uninformative prior
# X_fixed_effect_mean = matrix(rep(0, times = n_parameters * (n_treatments+1)), 
#                            nrow = n_parameters, ncol = n_treatments+1, byrow = T)
# X_fixed_effect_precision = diag(n_treatments+1)
# parameter_rate_matrix = diag(.1,n_parameters)

# ---- Effect of M on Y mean and precision ----
# immensely informative priors
M_fixed_effect_mean = t(mediator_effect_matrix)
M_fixed_effect_precision = diag(n_parameters)
Y_rate_matrix = solve(Y_covariance)

# # an uninformative prior
# M_fixed_effect_mean = matrix(rep(0, times = n_parameters * n_outcomes), 
#                              nrow = n_outcomes, ncol = n_parameters)
# M_fixed_effect_precision = diag(n_parameters)
# 
# Y_rate_matrix = diag(.1,n_outcomes)

# ---- Direct effect mean and precision priors ----
# immensely informative priors
direct_effect_mean = t(direct_effect)
direct_effect_precision = diag(n_treatments+1)

# # an uninformative prior
# direct_effect_mean = matrix(rep(0, times = n_treatments+1),
#                             nrow = n_outcomes, ncol = n_treatments+1)
# direct_effect_precision = diag(n_treatments+1)

#checking plots 
plot(Y, X[,2])
plot(Y, t(parameter_matrix[3,]))

jags_data = list(X = X,
                 parameter_matrix = parameter_matrix,
                 Y = Y,
                 nsubs = nsubs,
                 n_treatments = n_treatments+1,
                 n_parameters = n_parameters,
                 n_mediators = n_mediators,
                 n_outcomes = n_outcomes,
                 X_fixed_effect_mean = X_fixed_effect_mean,
                 X_fixed_effect_precision = X_fixed_effect_precision,
                 M_fixed_effect_mean = M_fixed_effect_mean,
                 M_fixed_effect_precision = M_fixed_effect_precision,
                 direct_effect_mean = direct_effect_mean,
                 direct_effect_precision = direct_effect_precision, 
                 parameter_rate_matrix = parameter_rate_matrix,
                 n_mediators_squared = n_mediators^2) # because JAGS doesn't recognize n_mediators^2 as an integer when used in a sequence operator

initial_values = list(X_fixed_effect = treatment_effect_matrix,
                      M_fixed_effect = mediator_effect_matrix,
                      indirect_effect = indirect_effect,
                      direct_effect = direct_effect)
# this should be adistribution
inits1 = list(.RNG.name = "base::Wichmann-Hill", 
              .RNG.seed = the_seed)
inits2 = list(.RNG.name = "base::Wichmann-Hill", 
              .RNG.seed = the_seed+500)

jagsModel = jags.model(file = "jags_model.R", data = jags_data, inits=list(inits1,inits2), n.chains = 2, n.adapt = 4000) #n.adapt = 4000
update(jagsModel, n.iter = 2500)

# warm up/burn in period 2500
# thinning

parameterlist = c("X_fixed_effect", "M_fixed_effect", "indirect_effect", "direct_effect")

codaSamples = coda.samples(jagsModel, variable.names = parameterlist, n.iter = 10000, thin = 1) 

resulttable = zcalc(codaSamples)
resulttable

plot(codaSamples[[1]])

source('postcalc.R')
zcalc(simple_chain)

# the dungeon
library(BayesianMediationA)
test = bma.bx.cy(pred = X[,-1],m = t(parameter_matrix),y =  Y, n.iter=1000,n.burnin = 1)
summary(test)

for(this_column in 1:n_mediators){
  for(this_row in 1:n_mediators){
    print((this_column-1)*n_mediators+this_row+n_mediators)
  }
}
