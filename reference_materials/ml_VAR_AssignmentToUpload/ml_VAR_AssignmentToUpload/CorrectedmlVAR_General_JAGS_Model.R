
model{
  
  # ---- (1) Likelihood ----
  # First, we write out the likelihood function
  for (pp in 1:P) {
    # 1st observation
    Y[pp,1,1:D] ~ dmnorm(mus1,prec1)
    
    # Recall that the ML-VAR(1) model looks like this (see Eq 3; Li et al., 2022)
    # y_{i,t} = Phi_{0,i} + Phi_{1,i}*(y_{i,t-1} - Phi_{0,i}) + u_{i,t}
    # You can also write as:
    # (y_{i,t} - Phi_{0,i}) = Phi_{0,i} + Phi_{1,i}*(y_{i,t-1} - Phi_{0,i}) + u_{i,t}
    # How fast I return to my person-specific baseline (or intercept), Phi_{0,i} at time t
    # depends on how far I was from my baseline, (y_{i,t-1} - Phi_{0,i}) at time t.
    
    
    # # ---- y under missingness ----
    # # This part uses the values implied by the fitted model to impute or "fill in" missing values
    # for(timePoint in Tmiss[pp,1:nmiss[pp]]){
    #   mus[pp,timePoint,1:D] <- intercept[pp,1:D] + AR[pp,1:D,1:D] %*% (Y[pp, timePoint-1,1:D] - intercept[pp,1:D])
    #   
    #   # NB: Here we use 2 univariate normal distributions for y1-yD
    #   for (dd in 1:D){
    #     Y[pp, timePoint, dd] ~ dnorm(mus[pp, timePoint, dd], pow(sd_noise[dd],-2))
    #   }
    # }
    
    # ---- y without missingness ----
    # This part pertains to data points that are NOT missing
    for(timePoint in Tseen[pp,2:nseen[pp]]){
      mus[pp,timePoint,1:D] <- intercept[pp,1:D] + AR[pp,1:D,1:D] %*% (Y[pp, timePoint-1,1:D] - intercept[pp,1:D])
      
      # NB: Here we use a multivariate normal distribution
      Y[pp,timePoint,1:D] ~ dmnorm(mus[pp,timePoint,1:D], prec_VAR[1:D,1:D]) 
    }
    
    
    # ----  Person-specific parameters (see Eq 6; Li et al., 2022) ----
    for (dd in 1:D){
      
      # mu1-mu2
      muInter[pp,dd] <- sum(Coefmus[1:nrCoeffmus,dd]*personCovariates[pp,1:nrCoeffmus])    
      
      # Intercepts of ARs and CRs - a1, a2, b1, b2
      for (dd2 in 1:D){ 
        muAR[pp,dd,dd2] <- sum(CoeffAR[1:nrCoeffAR,dd,dd2]*personCovariates[pp,1:nrCoeffAR])  
      } # close loop over dd2 (dimensions)
      
    } # close loop over dd (dimensions)
    
    # ---- Multivariate distribution for correlated person-specific parameters ----
    
    # Now assembling mu1, mu2, a1, & a2 into bmu to define the mean vector of a 
    # multivariate normal distribution (MVM) of the correlated person-specific parameters
    bmu[pp,1:D] = muInter[pp,1:D] #@@@ nrCoeffmus --> D
    for (dd in 1:D){
      bmu[pp,(D+dd)] = muAR[pp,dd,dd] #@@@ nrCoeffmus --> D
    }
    # Now we specify the person-specific mu1-a2 as multivariate normally distributed
    b[pp,1:nrW] ~ dmnorm(bmu[pp,1:nrW],bpre[1:nrW,1:nrW]) #bpre is the precision matrix for the correlated random effects
    
    #Gather person-specific parameters into matrices with more intuitive names for easier tracking later
    #mu_1i - mu_2i
    for (dd in 1:D){
      intercept[pp,dd] = b[pp,dd]
      
      # AR
      AR[pp,dd,dd] = b[pp, D+dd]#@@@ nrCoeffmus --> D
      
    } #Close dd loop
    
    # ---- Univariate dist for uncorrelated person-specific parameters ----
    for (dd in 1:(D-1)) {
      for (dd2 in (dd+1):D) {
        # CR - univariate normal distributions to allow the random effects for the 
        # CR coefs to be uncorrelated with other person-specific parameters
        AR[pp, dd, dd2] ~ dnorm(muAR[pp,dd,dd2],pow(sigmaAR[dd, dd2],-2)) #Person-specific CR y2 -> y1
        AR[pp, dd2, dd] ~ dnorm(muAR[pp,dd2,dd],pow(sigmaAR[dd2, dd],-2)) #Person-specific CR y1 -> y2
      } #Close dd2 loop
    }#dd loop
    
  } #Close loop over persons
  
  # ---- Person-invariant process noise cov structure ----
  
  # lsigma1, lsigma2
  for (dd in 1:D){
    logS_Means[dd] <- CoefflogSD_VAR[dd]   
    
    # Exponentiating the log so we get the person-specific process noise SD for each DV
    sd_noise[dd] <- exp(logS_Means[dd])
    
    # Covariance matrix, Sigma_VAR, for the process noises. Notice it uses the correlation for the process
    # noise (corr_noise[pp]) here. The computations below basically allow you to take the univariate,
    # person-specific SDs and correlations of the process noises, and create the process noise cov matrix
    Sigma_VAR[dd,dd] = sd_noise[dd] * sd_noise[dd]
    
    # Fisher transform mean of correlation between process noises (z) = Coeffcorr_VAR
    # Perform transformation from z back to r to get the person-specific correlations between process noises
    corr_ij[dd] <- (exp(2*Coeffcorr_VAR[dd])-1)/ (exp(2*Coeffcorr_VAR[dd])+1) 
  } #Close dd loop
  
  
  for (dd in 1:(D-1)) {
    for (dd2 in (dd+1):D) {
      v <- sd_noise[dd] * corr_ij[dd] * sd_noise[dd2] #to check
      Sigma_VAR[dd, dd2] <- v
      Sigma_VAR[dd2, dd] <- v
    }
  }
  
  # Do remember to invert the process noise cov matrix to give jags the precision matrix
  prec_VAR[1:D,1:D] <- inverse(Sigma_VAR[1:D,1:D])
  
  
  # ---- (2) Priors ----
  
  # ---- (2a) Fixed effects coefs for the person-specific parameters ----
  for (cov in 1:nrCoeffmus){
    for (dd in 1:D){
      
      # mu1-mu2
      Coefmus[cov,dd] ~ dnorm(0, 1) #THINK: Is a standard normal distribution a good prior for your data set points?
    }
  }
  
  for (cov in 1:nrCoeffAR){
    for (dd in 1:D){
      for (dd2 in 1:D) {
        
        # Fixed effects coefficients for a1-b2
        CoeffAR[cov,dd,dd2] ~ dnorm(0, 1) #THINK: This may be quite diffuse for the AR and CR coefs, but depends on scales of your predictors
      }
    }
  }
  
  # ---- (2b) Process noise-related priors ----
  
  # lsigma1-lsigma2, z_r of process noises
  for (dd in 1:D){
    CoefflogSD_VAR[dd] ~ dnorm(0, 1) #THINK: Standard normal for log of SD(process noise) -> (exp(z))^2 in terms of process noise variances
    
    #z-transformed correlation between process noises
    Coeffcorr_VAR[dd] ~ dnorm(0, 1) #Standard normal - should be quite decent for z scores
  }
  
  
  # ---- (2c) Stuff related to person-specific parameters w/ uncorrelated random effects (b1-b2,lsigma1-2) ----
  # Standard deviations of random effects (see Psi in Eq 6; Li et al., 2022)
  for (dd in 1:D) {
    logS_SD[dd] ~ dunif(0, 1) #Log transformations of SDs of random effects for lsigma1-lsigma2 -> Psi[dd] = (exp(logS_SD))^2
  }
  for (dd in 1:(D-1)) {
    for (dd2 in (dd+1):D) {
      sigmaAR[dd,dd2] ~ dunif(0, 1) # Random effect SD for cross-regression from y_dd2 to y_dd
      sigmaAR[dd2,dd] ~ dunif(0, 1) # Random effect SD for cross-regression from y_dd to y_dd2
    }
  }
  
  # ---- (2d) Stuff related to person-specific parameters w/ correlated random effects (mu1-a_2) ----
  # Prior for inverse of the random effect covariance matrix as following a Wishart distribution
  # In other words, the random effect covariance matrix follows an inverse-Wishart (IW) distribution
  # In a Wishart distribution, you have to specify the a scale matrix (basically a prior mean covariance 
  # you believe to be plausible) and degrees of freedom (strength of your belief) - small values like 1-5
  # are weakly informative
  # Li et al claimed they used an identity matrix as the scale matrix (W), and df of number of random effects + 3
  # But in the jags main set-up code (MLVAR_JAGS_Code.R), you can see that they set the scaling matrix as
  # the scale matrix of the IW distribution: W = diag(c(0.5,0.5,0.1,0.1),4)
  
  bpre[1:nrW,1:nrW] ~ dwish(W,nrW+3) # Precision matrix for random effects follow a Wishart
  
  bcov[1:nrW,1:nrW] <- inverse(bpre[1:nrW,1:nrW]) # Get random effects covariance matrix by taking the inverse
  
  # Now store the SDs of the random effects for mu1-a2.
  for (dd in 1:D){
    sigmaInter[dd] = sqrt(bcov[dd,dd])
    sigmaAR[dd,dd] = sqrt(bcov[D+dd,nrCoeffmus+dd])#@@@ nrCoeffmus --> D
  }
}# End of model

