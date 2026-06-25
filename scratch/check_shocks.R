library(dplyr)
library(purrr)
library(readr)

devtools::load_all(here::here())

# 1. Load transition counts and patient risk profiles
data("ipd_risk")

# 2. Load the BUGS MCMC output
load(here::here("data/jags_output.RData"))

n_total <- nrow(ipd_risk)
interv_names <- c("obs", "risk6", "risk4")
n_sim <- 1 # For deterministic base-case comparison, we only need 1 simulation
J <- 12 # 12-year time horizon

# Setup initial starting populations from individual patient data
init <- inits_by_interv(ipd_risk, group_names = c("icd", "risk_over_6", "risk_over_4"))
init_counts <- list(
  obs = list(
    ICD    = dplyr::filter(init, risk == "icd", status == "ICD")$n,
    no_ICD = dplyr::filter(init, risk == "icd", status == "low_risk")$n
  ),
  risk6 = list(
    ICD    = dplyr::filter(init, risk == "risk_over_6", status == "ICD")$n,
    no_ICD = dplyr::filter(init, risk == "risk_over_6", status == "low_risk")$n
  ),
  risk4 = list(
    ICD    = dplyr::filter(init, risk == "risk_over_4", status == "ICD")$n,
    no_ICD = dplyr::filter(init, risk == "risk_over_4", status == "low_risk")$n
  )
)

S <- 6
n_init <- list(
  obs    = c(init_counts$obs$ICD,    0, 0, init_counts$obs$no_ICD,    0, 0),
  risk6  = c(init_counts$risk6$ICD,  0, 0, init_counts$risk6$no_ICD,  0, 0),
  risk4  = c(init_counts$risk4$ICD,  0, 0, init_counts$risk4$no_ICD,  0, 0)
)

STATE_ICD   <- 1
STATE_SHOCK <- 2
STATE_DEATH <- 3

# Compute mean transition matrix from BUGS output
lambda.0_mean <- apply(mm1$BUGSoutput$sims.list$lambda.0, c(2, 3), mean)
lambda.1_mean <- apply(mm1$BUGSoutput$sims.list$lambda.1, c(2, 3), mean)
lambda.2_mean <- apply(mm1$BUGSoutput$sims.list$lambda.2, c(2, 3), mean)

# Copy background mortality (ICD -> Death) to Shock state (Shock -> Death)
lambda.0_mean[STATE_SHOCK, STATE_DEATH] <- lambda.0_mean[STATE_ICD, STATE_DEATH]
lambda.1_mean[STATE_SHOCK, STATE_DEATH] <- lambda.1_mean[STATE_ICD, STATE_DEATH]
lambda.2_mean[STATE_SHOCK, STATE_DEATH] <- lambda.2_mean[STATE_ICD, STATE_DEATH]

# Set Shock -> ICD transition as the remainder
lambda.0_mean[STATE_SHOCK, STATE_ICD]   <- 1 - lambda.0_mean[STATE_SHOCK, STATE_DEATH]
lambda.1_mean[STATE_SHOCK, STATE_ICD]   <- 1 - lambda.1_mean[STATE_SHOCK, STATE_DEATH]
lambda.2_mean[STATE_SHOCK, STATE_ICD]   <- 1 - lambda.2_mean[STATE_SHOCK, STATE_DEATH]

# Probability transition lists
probs_empty <- array(NA, dim = c(S, S, 1, n_sim))
probs <- list(obs = probs_empty, risk6 = probs_empty, risk4 = probs_empty)

probs$obs[,, 1, 1]    <- lambda.0_mean
probs$risk6[,, 1, 1]  <- lambda.1_mean
probs$risk4[,, 1, 1]  <- lambda.2_mean

pop <- list()
for (i in interv_names) {
  pop[[i]] <- array(NA, dim = c(S, J, n_sim))
  pop[[i]][, 1, 1] <- n_init[[i]]
}

# Run simulation
for (k in interv_names) {
  for (j in 2:J) {
    for (s in 1:S) {
      pop[[k]][s, j, 1] <- t(pop[[k]][, j-1, 1]) %*% probs[[k]][, s, 1, 1]
    }
  }
}

# Cumulative shocks is the sum of people entering the Shock state (State 2) from year 2 to 12,
# plus those who start in the Shock state (which is 0).
# Let's print the sum of populations in the Shock state at each cycle.
# Note that in a discrete Markov model, the population in State 2 at cycle j is the number of shocks in cycle j.
shocks_obs <- sum(pop$obs[2, , 1])
shocks_risk6 <- sum(pop$risk6[2, , 1])
shocks_risk4 <- sum(pop$risk4[2, , 1])

cat("Observed Strategy cumulative shocks:", shocks_obs, "\n")
cat("risk6 Strategy cumulative shocks:", shocks_risk6, "\n")
cat("risk4 Strategy cumulative shocks:", shocks_risk4, "\n")
