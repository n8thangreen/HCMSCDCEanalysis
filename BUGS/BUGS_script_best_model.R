# run BUGS model script
# Combined Best Model: Continuous Poisson Rates with Bivariate Normal prior on correlated death rates
#
# estimate transition probabilities
# from discrete time counts
# using continuous Poisson rates and converting to probabilities

library(R2jags)
library(here)
# devtools::load_all(".") # load package for transdat data

# Load aggregated transition counts from the package
# created by scripts/prep_study_data.R
data("transdat")

data_obs <- transdat$icd
data_risk6 <- transdat$risk_over_6
data_risk4 <- transdat$risk_over_4
data_no_icd <- transdat$no_icd

n_S <- nrow(data_obs$ICD)  # number of states

# Pad with 0s to align matrices
empty_mat <- matrix(0, nrow = n_S, ncol = n_S)

r.0 <- rbind(cbind(data_obs$ICD[, 1:n_S], empty_mat),
             cbind(empty_mat, data_obs$low_risk[, 1:n_S]))

r.1 <- rbind(cbind(data_risk6$ICD[, 1:n_S], empty_mat),
             cbind(empty_mat, data_risk6$low_risk[, 1:n_S]))

r.2 <- rbind(cbind(data_risk4$ICD[, 1:n_S], empty_mat),
             cbind(empty_mat, data_risk4$low_risk[, 1:n_S]))

r.3 <- rbind(cbind(data_no_icd$ICD[, 1:n_S], empty_mat),
             cbind(empty_mat, data_no_icd$low_risk[, 1:n_S]))

n.0 <- unname(rowSums(r.0))
n.1 <- unname(rowSums(r.1))
n.2 <- unname(rowSums(r.2))
n.3 <- unname(rowSums(r.3))

# Combine strategy inputs into vectors for JAGS
# 1 = observed, 2 = >6% risk, 3 = >4% risk, 4 = no ICD
n_icd <- c(n.0[1], n.1[1], n.2[1], n.3[1])
n_noicd <- c(n.0[4], n.1[4], n.2[4], n.3[4])

r_icd_shock <- c(r.0[1, 2], r.1[1, 2], r.2[1, 2], r.3[1, 2])
r_icd_death <- c(r.0[1, 3], r.1[1, 3], r.2[1, 3], r.3[1, 3])
r_noicd_scd <- c(r.0[4, 5], r.1[4, 5], r.2[4, 5], r.3[4, 5])
r_noicd_death <- c(r.0[4, 6], r.1[4, 6], r.2[4, 6], r.3[4, 6])

# Specify prior values for Bivariate Normal (correlated death rates)
mean_prior <- c(-5, -5)              # Mean log-rates of death (roughly 0.0067)
prec_prior <- diag(2) * 0.1          # Uninformative precision matrix
R_wishart <- diag(2) * 0.1           # Scale matrix for Wishart prior on precision

dataJags <- list(
  n_icd = n_icd,
  n_noicd = n_noicd,
  r_icd_shock = r_icd_shock,
  r_icd_death = r_icd_death,
  r_noicd_scd = r_noicd_scd,
  r_noicd_death = r_noicd_death,
  mean_prior = mean_prior,
  prec_prior = prec_prior,
  R_wishart = R_wishart,

  #
  from_shock = c(1, 0, 0, 0, 0, 0),
  from_ICD_death = c(0, 0, 1, 0, 0, 0)
)

filein <- "BUGS/BUGS_code_best_model.txt"

params <- c("lambda.0", "lambda.1", "lambda.2", "lambda.3",
            "rate_icd_shock", "rate_icd_death", "rate_noicd_scd", "rate_noicd_death")

# Initialize rates for starting values
inits <- function() {
  rate_init <- array(NA, dim = c(4, 4, 6))
  for (i in 1:4) {
    rate_init[i, 1, 2] <- 0.05
    rate_init[i, 4, 5] <- 0.01
  }
  list(
    rate = rate_init,
    log_rates_death = matrix(rnorm(8, -5, 0.5), nrow = 4, ncol = 2)
  )
}

n.iter <- 10000
n.burnin <- 5000
n.thin <- floor((n.iter - n.burnin) / 500)

mm1 <- jags(
  data = dataJags,
  inits = inits,
  parameters.to.save = params,
  model.file = filein,
  n.chains = 2,
  n.iter = n.iter,
  n.burnin = n.burnin,
  n.thin = n.thin,
  DIC = TRUE
)

# Attach MCMC simulations and export
R2WinBUGS::attach.bugs(mm1$BUGSoutput)

# Save output files to be consumed by the CE Markov simulation scripts
save(mm1, file = here::here("data/jags_output.RData"))
saveRDS(lambda.0, file = here::here("data/lambda0.Rds"))
saveRDS(lambda.1, file = here::here("data/lambda1.Rds"))
saveRDS(lambda.2, file = here::here("data/lambda2.Rds"))
saveRDS(lambda.3, file = here::here("data/lambda3.Rds"))
