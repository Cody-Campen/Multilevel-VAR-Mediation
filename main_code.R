#' This is the code which runs the JAGS model and calculates the answers
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

the_seed = the_seed

answers = NULL # this is the dataframe we will store the answers in

# Dataset generation variables
n_people = 138
n_times = 56
percent_missing = .2
n_treatments = 2
n_mediators = 2
n_parameters = 4
n_outcomes = 1
treatment_effect_matrix = matrix(c(0, 0, 1,   # X -> M_1 intercept
                                   0, 0, 0,   # X -> M_2 intercept
                                   0, .2, 0,   # X -> M_1 autoregression
                                   0, 0, 0,   # X -> M_1 to M_2 crossregression
                                   0, 0, 0,   # X -> M_2 to M_1 crossregression #-.2
                                   0, 0, 0),  # X -> M_2 autoregression
                                 nrow = 6, ncol = n_treatments+1, byrow = T)
mediator_effect_matrix = matrix(c(1,   # M_1 intercept -> Y
                                  0,   # M_2 intercept -> Y
                                  1,   # M_1 autoregression -> Y
                                  0,   # M_1 -> M_2 crossregression -> Y
                                  0,   # M_2 -> M_1 crossregression -> Y # 1
                                  0),  # M_2 autoregression -> Y
                                nrow = 6, ncol = n_outcomes, byrow = T)
direct_effect = matrix(c(0,  # Y intercept
                         1,  # X1 -> Y
                         -1), # X2 -> Y
                       nrow = n_treatments + 1, ncol = n_outcomes, byrow = T)
parameter_matrix_covariance = diag(c(10, 10, .1, .1, .1, .1))
Y_covariance = diag(n_outcomes)/1
m_standard_deviations = c(1, 1)
m_correlations = matrix(c(1, .2,
                          .2, 1),
                        nrow = 2, ncol = 2)
m_covariance = diag(m_standard_deviations) %*% m_correlations %*% diag(m_standard_deviations)

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
M = array(as.matrix(dataset[,c("M1", "M2")]), dim = c(n_times, n_people, n_mediators)) |> aperm(c(2,1,3))
# you want the dimensions to be structured: n_people, n_times, n_mediators

# Structure the outcome
Y = as.matrix(dataset[dataset$time == 1, c("Y1")], nrows = n_people)

# Now, we create objects to handle the time indexing
# Begin by creating objects to store our values
n_seen = vector(length = n_people)
n_miss = vector(length = n_people)
times_seen = matrix(NA, nrow = n_people, ncol = n_times)
times_missed = matrix(NA, nrow = n_people, ncol = n_times)

# and loop over each person to fill in their values based on the missingness index, is_missing.
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

# ---- Set up the initial values ----
inits_list = create_inits(n_chains = n_chains, 
                          n_treatments = n_treatments+1, 
                          n_mediators = n_mediators,
                          n_outcomes = n_outcomes)

# ---- Define the JAGS data ----
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
                 X_fixed_effect_mean = rep(0, times = n_treatments+1),
                 M_fixed_effect_mean = rep(0, times = n_parameters),
                 direct_effect_mean = rep(0, times = n_treatments+1),
                 X_to_intercept.precision = matrix(c(.01, 0, 0,
                                                     0, .01, 0,
                                                     0, 0, .01),
                                                   nrow = 3),
                 X_to_AR.precision = matrix(c(1, 0, 0,
                                              0, 1, 0,
                                              0, 0, 1),
                                            nrow = 3),
                 X_to_CR.precision = matrix(c(1, 0, 0,
                                              0, 1, 0,
                                              0, 0, 1),
                                            nrow = 3),
                 M_fixed_effect.precision = matrix(c(.1, 0, 0,
                                                     0, .1, 0,
                                                     0, 0, .1),
                                                   nrow = 3),
                 direct_effect.precision = matrix(c(.1, 0, 0,
                                                    0, .1, 0,
                                                    0, 0, .1),
                                                  nrow = 3),
                 parameter_matrix.rate = matrix(c(15, 0, 0, 0, 0, 0,
                                                  0, 15, 0, 0, 0, 0,
                                                  0, 0, .25, 0, 0, 0,
                                                  0, 0, 0, .25, 0, 0,
                                                  0, 0, 0, 0, .25, 0,
                                                  0, 0, 0, 0, 0, .25),
                                                nrow = 6))

# ---- Compile and run the model ---
jagsModel = jags.model(file = "jags_model.R", 
                       data = jags_data, 
                       # inits = inits_list, 
                       n.chains = n_chains, 
                       n.adapt = 4000)

update(jagsModel, n.iter = 6000)

parameterlist = c("X_to_intercept",
                  "X_to_AR",
                  "X_to_CR",
                  "parameter_matrix.precision",
                  "intercept_to_Y",
                  "AR_to_Y",
                  "CR_to_Y",
                  "direct_effect",
                  "Y.precision",
                  "intercept_indirect_effect",
                  "AR_indirect_effect",
                  "CR_indirect_effect",
                  "M.precision")

before_time = Sys.time()
codaSamples = coda.samples(jagsModel, variable.names = parameterlist, n.iter = 8000, thin = 1) 
run_time = Sys.time() - before_time

sample_names = varnames(codaSamples)

ordered_parameters = c()
for(this_parameter in parameterlist){
  matching_parameters = grep(paste0("^", this_parameter, "(\\[|$)"), sample_names, value = TRUE)
  ordered_parameters = c(ordered_parameters, matching_parameters)
}

# ---- Collect our simulation performance statistics ----

resulttable = zcalc(codaSamples[, ordered_parameters])

# we're going to do this hard coded for now, just so we can submit this job before the end of day.
true_values = c(treatment_effect_matrix[1:2, ],
                treatment_effect_matrix[c(3,6), ],
                # treatment_effect_matrix[4:5, ],
                parameter_matrix_covariance[1:n_parameters, 1:n_parameters]^(-2),
                mediator_effect_matrix[1:2, ],
                mediator_effect_matrix[c(3, 6), ],
                # mediator_effect_matrix[4:5, ],
                direct_effect,
                Y_covariance^(-2),
                treatment_effect_matrix[1:2, ] * mediator_effect_matrix[1:2,],
                treatment_effect_matrix[3:4, ] * mediator_effect_matrix[c(3, 6), ],
                # treatment_effect_matrix[5:6, ] * mediator_effect_matrix[4:5, ],
                m_covariance^(-2))
                
raw_bias = resulttable$mean - true_values
names(raw_bias) = paste(rownames(resulttable),"raw_bias", sep = "_")

relative_bias = (resulttable$mean - true_values) / (true_values)
names(relative_bias) = paste(rownames(resulttable),"relative_bias", sep = "_")

rmse = sqrt((resulttable$mean - true_values)^2)
names(rmse) = paste(rownames(resulttable), "rmse", sep = "_")

sd = resulttable$sd
names(sd) = paste(rownames(resulttable), "sd", sep = "_")

CI_low = resulttable$`CrI 2.5%`
names(CI_low) = paste(rownames(resulttable), "CI_low", sep = "_")

CI_high = resulttable$`CrI 97.5%`
names(CI_high) = paste(rownames(resulttable), "CI_high", sep = "_")

coverage = (resulttable$`CrI 2.5%` < true_values & true_values < resulttable$`CrI 97.5%`)
names(coverage) = paste(rownames(resulttable), "coverage", sep = "_")

HDI_low = resulttable$`95% HDI_L`
names(HDI_low) = paste(rownames(resulttable), "HDI_low", sep = "_")

HDI_high = resulttable$`95% HDI_H`
names(HDI_high) = paste(rownames(resulttable), "HDI_high", sep = "_")

ess = resulttable$ESS
names(ess) = paste(rownames(resulttable), "ess", sep = "_")

rhat = resulttable$RHAT
names(rhat) = paste(rownames(resulttable), "rhat", sep = "_")

answers = rbind(answers, c(this_sim = the_seed,
                           run_time = run_time,
                           raw_bias,
                           relative_bias,
                           rmse,
                           sd,
                           CI_low,
                           CI_high,
                           coverage,
                           HDI_low,
                           HDI_high,
                           ess,
                           rhat))

save(answers, file = paste0("sim_", the_seed, ".RData"))
