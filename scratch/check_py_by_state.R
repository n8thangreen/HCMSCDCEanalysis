library(dplyr)
library(purrr)
library(readr)

devtools::load_all(here::here())

# 1. Load transition counts and patient risk profiles
load(here::here("data/transdat.RData"))
data("ipd_risk")

# 2. Load the BUGS MCMC output
load(here::here("data/jags_output_all_strat_direct.RData"))
lambda_samples <- mm1$BUGSoutput$sims.list$lambda

# Get the mean transitions
lambda_mean <- apply(lambda_samples, c(2, 3, 4), mean)

n_total <- nrow(ipd_risk)
interv_names <- names(transdat)
n_sim <- 1
J <- 12 # 12-year time horizon
S <- 3  # state dimension in BUGS (hcm, scd/shock, death)

# Process probabilities to 2*S = 6 states total
probs <- list()
p_death_bg <- lambda_mean[, 1, 3] # ICD -> Death
lambda_mean[, 2, 3] <- p_death_bg
lambda_mean[, 2, 1] <- 1 - p_death_bg

for (j in seq_along(interv_names)) {
  nm <- interv_names[j]
  probs_array <- array(NA, dim = c(2 * S, 2 * S, 1, n_sim))
  probs_array[,, 1, 1] <- lambda_mean[j, , ]
  probs[[nm]] <- probs_array
}

# Setup initial populations
init_counts <- inits_by_interv(ipd_risk, group_names = interv_names) |>
  mutate(state = "hcm")

init_full <- expand.grid(
  risk = interv_names,
  status = c("ICD", "low_risk"),
  state = c("hcm", "scd/shock", "death")) |>
  arrange(risk, status) |>
  left_join(init_counts, by = c("risk", "status", "state")) |>
  tidyr::replace_na(list(n = 0))

pop0 <- split(init_full$n, init_full$risk)
pop <- list()
for (i in interv_names) {
  pop[[i]] <- array(NA, dim = c(2*S, J, n_sim))
  pop[[i]][, 1, 1] <- pop0[[i]]
}

# Run simulation
for (k in interv_names) {
  for (j in 2:J) {
    for (s in 1:(2*S)) {
      pop[[k]][s, j, 1] <- t(pop[[k]][, j-1, 1]) %*% probs[[k]][, s, 1, 1]
    }
  }
}

# Discount factor (3.5% per year)
disc <- 1.035^(0:11)

for (nm in c("risk_over_4", "risk_over_6", "icd")) {
  cat("\n=== STRATEGY:", nm, "===\n")
  # Discounted person-years in each state
  icd_py <- sum(pop[[nm]][1, , 1] / disc) / n_total
  shock_py <- sum(pop[[nm]][2, , 1] / disc) / n_total
  noicd_py <- sum(pop[[nm]][4, , 1] / disc) / n_total
  scd_py <- sum(pop[[nm]][5, , 1] / disc) / n_total
  
  cat("Discounted PY in ICD state (1):", icd_py, "\n")
  cat("Discounted PY in Shock state (2):", shock_py, "\n")
  cat("Discounted PY in No-ICD state (4):", noicd_py, "\n")
  cat("Discounted PY in SCD state (5):", scd_py, "\n")
}
