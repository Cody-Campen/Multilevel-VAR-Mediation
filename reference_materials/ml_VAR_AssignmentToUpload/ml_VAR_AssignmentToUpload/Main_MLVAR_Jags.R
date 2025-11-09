rm(list=ls())
library(rjags)
library(dplyr)
source("postcalc.R")

r1=1  #Factor to set random seed number
N=100 #Number of participants
O=60  #Number of time points per person
D=2 #Number of DVs
nrCovariates = 2 #Number of person-specific covariates including "1"

# ---- Read in data ----
VAR_DemoData = read.table("SimulatedData_low_N100T60_InvSigma.dat",na.strings="99")
colnames(VAR_DemoData) = c("id","x","y1","y2")

# ---- Structure Y_obs and personCovariates ----
# Get Y_Obs data into dim (O,N,D), then permute so it becomes (N,O,D)
Y_obs <- array(as.matrix(VAR_DemoData[,c("y1","y2")]), dim = c(O, N, D)) %>% aperm(c(2,1,3))
dim(Y_obs)  # (N, O, D)

#Person-specific covariates
Time = rep(1:O,N)
X = cbind(rep(1,N), as.matrix(VAR_DemoData[Time==1,"x"], nrows=N))

# ---- Create missingness indicators ----
Index = matrix(NA,nrow=N,ncol=O)
for(pp in 1:N){
  for (oo in 1:O){
    Index[pp,oo] = ifelse(sum(is.na(Y_obs[pp,oo,]))==0, 0, 1)
  }
}

#Create a matrix of NAs, with non-NA entries
Tmiss = matrix(NA,nrow=N,ncol=O)
Tseen = matrix(NA,nrow=N,ncol=O)
nmiss=c()
nseen=c()
for(pp in 1:N){
  nmiss[pp] = length(which(Index[pp,] %in% 1))
  if(nmiss[pp] != 0){
    Tmiss[pp, 1:nmiss[pp]] <- which(Index[pp,] %in% 1) #location of records that have at least one missing entry
  }
  nseen[pp] = length(which(Index[pp,] %in% 0))
  Tseen[pp, 1:nseen[pp]] <- which(Index[pp,] %in% 0)  #location of fully observed data
}

# ---- Create other specifications needed for priors ----
# the scale matrix of the IW distribution - # of elements = # of correlated person-specific parameters
W = diag(c(0.5,0.5,0.1,0.1),4)
# prior for the 1st observation
mus1 = c(0,0) #Means for y1(t=1) and y2(t=1)
prec1 = diag(2) #Precision (1/var^2) for y1(t=1) and y2(t=1)



############################ JAGS Codes################################################

jags_data <- list(Y = Y_obs, #Observed Y data
             personCovariates = X, #Matrix with person-specific covariates, including a column of "1s" 
             nrCoeffmus= dim(X)[2],  #Number of covariates to explain person-specific mus
             nrCoeffAR = dim(X)[2], #Number of covariates to explain person-specific AR
             P = dim(Y_obs)[1], #Number of participants
             D = dim(Y_obs)[3], #Number of DVs
             nmiss = nmiss, #Vector of number of missing occasions by participant
             nseen = nseen, #Vector of number of complete data occasions by participant
             Tmiss = Tmiss, #Matrix of location for each instance of missing data for each participant
             Tseen = Tseen, #Matrix of location for each instance of complete data for each participant
             W = W, mus1=mus1, prec1=prec1,
             nrW = dim(W)[1]
             )

inits1 <- list(.RNG.name="base::Wichmann-Hill", .RNG.seed=r1)
inits2 <- list(.RNG.name="base::Wichmann-Hill", .RNG.seed=r1+500)


jagsModel <- jags.model(file = "CorrectedmlVAR_General_JAGS_Model.R", data = jags_data, inits=list(inits1,inits2), n.chains = 2, n.adapt = 1000) #n.adapt = 4000
update(jagsModel, n.iter = 1000) 


parameterlist <-c("CoeffInter","CoeffAR","CoefflogSD_VAR","Coeffcorr_VAR",
                  "sigmaInter","sigmaAR","logS_SD","R_SD","bcov",
                  "intercept", "AR", "corr_noise", "Sigma_VAR")
codaSamples <- coda.samples(jagsModel, variable.names = parameterlist, n.iter = 1000, thin = 1) 


resulttable <- zcalc(codaSamples)

#What values of a_hat and b_hat did you get? How do they compare to the true values?
plot(codaSamples)
gelman.plot(codaSamples)

#Fixed effects for AR and CR
resulttable[grep("CoeffAR",row.names(resulttable)),] 

#Fixed effects for mu1-mu2
resulttable[grep("CoeffInter",row.names(resulttable)),] 

warns = warnings()

## save results

save(warns, resulttable, file=paste0("InvSigma_Result_JAGS_low_N",N,"T",O,".Rdata"))




