# Example from John M. White's website at:
# https://www.johnmyleswhite.com/notebook/2010/08/20/using-jags-in-r-with-the-rjags-package/
# Last modified: 8/25/25
# Authors: Sy-Miin Chow & Jek Pesigan


rm(list = ls()) #Remove everything in your global environment - make sure you want to do that

# ---- Load library ----
library('rjags') #Needed for some of the diagnostic functions
library('R2jags') #Need for parallel execution of model fitting
#You can also load R2jags and coda without loading rjags
source('postcalc.R')


# ---- Example 1: Generate linear regression data ----
set.seed(123456)
N <- 1000
x <- 1:N
epsilon <- rnorm(N, 0, 1)
a = 0
b = 1
y <- a + b*x + epsilon

write.table(
  data.frame(X = x, Y = y, Epsilon = epsilon),
  file = 'example1.data',
  row.names = FALSE,
  col.names = TRUE
)

# ---- Example 1: Compile and run jags model ----
jags <- jags.model(
  'Model01_LinearRegression.txt',
  data = list(
    'x' = x,
    'y' = y,
    'N' = N
  ),
  #inits = function () {#Optional function to provide initial values for parameters across chains
  #  list(
  #    'a' = rnorm(1, 0, 10),
  #    'b' = runif(1, .0001, 10)
  #  )
  #}
  n.chains = 4, #Number of independent chains
  n.adapt = 100 #Number of internal JAGS adaptation steps
)

update(jags, 1000) #Number of burn-in iterations to discard

#We could use jags.sample or coda.samples to sample from the posterior distribution.
#It is easier to call other diagnostic functions after coda.samples so we are going to use that.
#jags.samples(jags,
#             c('a', 'b'),
#            1000)


samples2 = coda.samples(jags,
                        c('a','b','sigma'),
                        1000)
#What values of a_hat and b_hat did you get? How do they compare to the true values?
plot(samples2)
gelman.plot(samples2)
gelman.diag(samples2)
zcalc(samples2)

# Alternatively, we can use jags.parallel to run the entire process in parallel utilizing
# parallel cores on your machine. 
fit <- R2jags::jags.parallel(
  data = list(
    'x' = x,
    'y' = y,
    'N' = N
  ),
  parameters.to.save = c("a", "b", "sigma"),
  model.file = 'Model01_LinearRegression.txt',
  n.chains = 4, #Number of independent chains
  n.iter = 2000, #Number of total iterations per chain (including burn-in)
  n.burnin = 1000, #Number of burn-in iterations
  #n.adapt = 100, #Cannot set the number of internal JAGS adaptation steps - default to 100
  n.thin = 1 #Thinning rate
)

print(fit)
# Convert to mcmc.list
samps <- as.mcmc(fit)   # First convert fit to mcmc.list


# (optional) keep only the parameters you saved
samps <- samps[, c("a", "b", "sigma")]

plot(samps)
gelman.plot(samps)
gelman.diag(samps)
zcalc(samps)


# ---- Example 2: Generate data from logistic regression model ----
set.seed(123456)
N <- 1000
x <- runif(N,-4,4)#runif(N,-40,40)
z <- 0.01 * x - 5
y <- sapply(1 / (1 + exp(-z)), function(p) {rbinom(1, 1, p)})

write.table(
  data.frame(X = x, Z = z, Y = y),
  file = 'example2.data',
  row.names = FALSE,
  col.names = TRUE
)


# ---- Example 2: Compile and run jags model ----
jags <- jags.model(
  'Model02_LogisticRegression.txt',
  data = list(
    'x' = x,
    'y' = y,
    'N' = N
  ),
  n.chains = 4,
  n.adapt = 100,
  inits = function () {#Optional function to provide initial values for parameters across chains
    list(
      'a' = rnorm(1, mean=0, sd=10),
      'b' = rnorm(1, mean=0, sd=10)
    )
  }
)

update(jags, 1000)

samples3 = coda.samples(
  jags,
  c('a', 'b'),
  1000
)

plot(samples3)
gelman.diag(samples3)

#DIY: Now let's try changing the data generation for x to runif(N,-40,40). What
#did you observe?