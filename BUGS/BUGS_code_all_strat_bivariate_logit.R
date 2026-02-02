model {
  # Modeling considering multi-parameter correlations

  # Hyperparameters for the multivariate normal distribution
  mu[1:2] ~ dmnorm(mu0[1:2], Omega[1:2,1:2])

  # Inverse of the covariance matrix (precision matrix)
  Omega[1:2,1:2] ~ dwish(R[1:2,1:2], 2)
  Sigma[1:2,1:2] <- inverse(Omega[1:2,1:2])

  # Between-study variation parameters
  for(i in 1:2) {
    tau[i] ~ dgamma(0.001, 0.001)
  }

  # Study-specific effects (multivariate random effects)
  for(k in 1:n_studies) {
    theta[k, 1:2] ~ dmnorm(mu[1:2], Omega[1:2,1:2])

    # Apply between-study variation
    log_or_shock[k] ~ dnorm(theta[k,1], tau[1])
    log_or_mortality[k] ~ dnorm(theta[k,2], tau[2])

    # Likelihood function (meta-analysis data)
    y_shock[k] ~ dnorm(log_or_shock[k], prec_shock[k])
    y_mortality[k] ~ dnorm(log_or_mortality[k], prec_mortality[k])
  }

  # Predictive distribution
  theta_pred[1:2] ~ dmnorm(mu[1:2], Omega[1:2,1:2])

  # Use observational data from Green et al. (2024) for baseline probabilities
  logit_p0_shock ~ dnorm(-5.986, 4)    # Corresponds to 0.25% baseline shock rate
  logit_p0_mortality ~ dnorm(-4.671, 4) # Corresponds to 0.93% baseline mortality rate
  logit_p0_scd ~ dnorm(-5.183, 4)      # Corresponds to 0.56% baseline SCD rate

  # Apply covariate effects
  logit_p_shock <- logit_p0_shock + theta_pred[1] * 0
  logit_p_mortality <- logit_p0_mortality + theta_pred[2] * 0
  logit_p_scd <- logit_p0_scd + 0  # Simplified handling

  # Transition probability matrix
  lambda[1,1] <- 1 - (ilogit(logit_p_shock) + ilogit(logit_p_mortality))
  lambda[1,2] <- ilogit(logit_p_shock)
  lambda[1,3] <- 0
  lambda[1,4] <- 0
  lambda[1,5] <- ilogit(logit_p_mortality)

  # State 2: ICD Shock
  p_recover_mpes ~ dbeta(75, 25)
  lambda[2,1] <- p_recover_mpes
  lambda[2,2] <- 1 - (p_recover_mpes + ilogit(logit_p_mortality) * 1.5)
  lambda[2,3] <- 0
  lambda[2,4] <- 0
  lambda[2,5] <- ilogit(logit_p_mortality) * 1.5

  # State 3: HCM without ICD
  lambda[3,1] <- 0
  lambda[3,2] <- 0
  lambda[3,3] <- 1 - (ilogit(logit_p_scd) + ilogit(logit_p_mortality) * 0.8)
  lambda[3,4] <- ilogit(logit_p_scd)
  lambda[3,5] <- ilogit(logit_p_mortality) * 0.8

  # States 4 and 5: Absorbing states
  lambda[4,1] <- 0
  lambda[4,2] <- 0
  lambda[4,3] <- 0
  lambda[4,4] <- 1
  lambda[4,5] <- 0

  lambda[5,1] <- 0
  lambda[5,2] <- 0
  lambda[5,3] <- 0
  lambda[5,4] <- 0
  lambda[5,5] <- 1
}
