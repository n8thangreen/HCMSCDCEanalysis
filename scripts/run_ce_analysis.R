# Unified Cost-Effectiveness Markov Simulation Runner
#
# Use this script to run the analysis for different devices (ICD vs. SICD)
# and model structures (Standard vs. Tunnel state models).

library(dplyr)
library(purrr)
library(BCEA)
library(readr)
library(ggplot2)
devtools::load_all(here::here())

# Set your configs ---

device <- "ICD"       # Options: "ICD" (Transvenous) or "SICD" (Subcutaneous)
model_type <- "standard" # Options: "standard" (6-state model) or "tunnel" (with replacement cycles)
J <- 12  #30                  # Time horizon / cycle limit (years)

# 1. Load transition counts and patient risk profiles
data("transdat")
data("ipd_risk")

# 2. Load the BUGS MCMC output
load(here::here("data/jags_output.RData"))
R2WinBUGS::attach.bugs(mm1$BUGSoutput) # makes lambda.0, lambda.1, lambda.2, lambda.3 available

n_total <- nrow(ipd_risk)
risk_thresh <- c(0.04, 0.06)
interv_names <- c("risk4", "risk6", "obs", "no_icd")
n_sim <- mm1$BUGSoutput$n.sims

# 3. Read baseline scenario parameters
inp_tab <- readr::read_csv(
  here::here("inst/extdata/scenario_input_values.csv"),
  col_names = TRUE, show_col_types = FALSE
)

# Map device to the corresponding scenario label in the CSV
lbl <- if (tolower(device) == "sicd") "sicd" else "tvicd"
inp <- inp_tab %>%
  dplyr::filter(Scenario_Label == lbl) %>%
  dplyr::slice(1)

# Dimension index definitions for array structuring
DIM_SIM  <- 1
DIM_FROM <- 2
DIM_TO   <- 3

TIME_SLICE_HOMOGENEOUS <- 1  # homogeneous transition probabilities only need 1 slice
YEAR_START             <- 1  # cohort starts at year 1

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
  ),
  no_icd = list(
    ICD    = 0,
    no_ICD = n_total
  )
)

# Model Structuring (Standard vs. Tunnel) ---

if (model_type == "standard") {
  S <- 6  # 6-state model
  inp$t_repl <- 1 # Enforce no tunnel states for standard 6-state model

  # Starting populations
  n_init <- list(
    obs    = c(init_counts$obs$ICD,    0, 0, init_counts$obs$no_ICD,    0, 0),
    risk6  = c(init_counts$risk6$ICD,  0, 0, init_counts$risk6$no_ICD,  0, 0),
    risk4  = c(init_counts$risk4$ICD,  0, 0, init_counts$risk4$no_ICD,  0, 0),
    no_icd = c(init_counts$no_icd$ICD, 0, 0, init_counts$no_icd$no_ICD, 0, 0)
  )

  # State index definitions for readability
  STATE_ICD   <- 1
  STATE_SHOCK <- 2
  STATE_DEATH <- 3

  # Copy background mortality (ICD -> Death) to Shock state (Shock -> Death)
  lambda.0[, STATE_SHOCK, STATE_DEATH] <- lambda.0[, STATE_ICD, STATE_DEATH]
  lambda.1[, STATE_SHOCK, STATE_DEATH] <- lambda.1[, STATE_ICD, STATE_DEATH]
  lambda.2[, STATE_SHOCK, STATE_DEATH] <- lambda.2[, STATE_ICD, STATE_DEATH]
  lambda.3[, STATE_SHOCK, STATE_DEATH] <- lambda.3[, STATE_ICD, STATE_DEATH]

  # Set Shock -> ICD transition as the remainder
  lambda.0[, STATE_SHOCK, STATE_ICD]   <- 1 - lambda.0[, STATE_SHOCK, STATE_DEATH]
  lambda.1[, STATE_SHOCK, STATE_ICD]   <- 1 - lambda.1[, STATE_SHOCK, STATE_DEATH]
  lambda.2[, STATE_SHOCK, STATE_ICD]   <- 1 - lambda.2[, STATE_SHOCK, STATE_DEATH]
  lambda.3[, STATE_SHOCK, STATE_ICD]   <- 1 - lambda.3[, STATE_SHOCK, STATE_DEATH]

  # Probability transition lists
  probs_empty <- array(NA, dim = c(S, S, TIME_SLICE_HOMOGENEOUS, n_sim))
  probs <- list(obs = probs_empty, risk6 = probs_empty, risk4 = probs_empty, no_icd = probs_empty)

  probs$obs[,, TIME_SLICE_HOMOGENEOUS, ]    <- aperm(lambda.0, c(DIM_FROM, DIM_TO, DIM_SIM))
  probs$risk6[,, TIME_SLICE_HOMOGENEOUS, ]  <- aperm(lambda.1, c(DIM_FROM, DIM_TO, DIM_SIM))
  probs$risk4[,, TIME_SLICE_HOMOGENEOUS, ]  <- aperm(lambda.2, c(DIM_FROM, DIM_TO, DIM_SIM))
  probs$no_icd[,, TIME_SLICE_HOMOGENEOUS, ] <- aperm(lambda.3, c(DIM_FROM, DIM_TO, DIM_SIM))

  # Populate starting cohort state populations
  pop <- list()
  for (i in interv_names) {
    pop[[i]] <- array(NA, dim = c(S, J, n_sim))
    pop[[i]][, YEAR_START, ] <- t(n_init[[i]])
  }

  # Run the standard 6-state simulation
  out <- ce_sim(
    pop = pop,
    probs = probs,
    inp = inp,
    c_unit = rc_unit,
    e_unit = re_unit,
    c_init = rc_init,
    e_init = re_init,
    distn = FALSE,
    interv_names = interv_names
  )

} else if (model_type == "tunnel") {
  n_tunnel <- inp$t_repl
  S <- 2 * n_tunnel + 4

  probs_empty <- array(NA, dim = c(S, S, TIME_SLICE_HOMOGENEOUS, n_sim))
  probs <- list(
    obs = probs_empty,
    risk6 = probs_empty,
    risk4 = probs_empty,
    no_icd = probs_empty
  )

  alpha <- 0
  for (i in 1:n_sim) {
    probs$obs[,, TIME_SLICE_HOMOGENEOUS, i]    <- lambda_tunnel(lambda.0[i,,], n_tunnel = n_tunnel, alpha)
    probs$risk6[,, TIME_SLICE_HOMOGENEOUS, i]  <- lambda_tunnel(lambda.1[i,,], n_tunnel = n_tunnel, alpha)
    probs$risk4[,, TIME_SLICE_HOMOGENEOUS, i]  <- lambda_tunnel(lambda.2[i,,], n_tunnel = n_tunnel, alpha)
    probs$no_icd[,, TIME_SLICE_HOMOGENEOUS, i] <- lambda_tunnel(lambda.3[i,,], n_tunnel = n_tunnel, alpha)
  }

  # Populate starting cohort state populations dynamically from data
  pop <- list()
  for (i in interv_names) {
    pop[[i]] <- array(NA, dim = c(S, J, n_sim))

    # Place initial patients in the first tunnel state (or healthy no-ICD state)
    pop[[i]][, YEAR_START, ] <- 0
    pop[[i]][1, YEAR_START, ] <- init_counts[[i]]$ICD
    pop[[i]][2 * n_tunnel + 2, YEAR_START, ] <- init_counts[[i]]$no_ICD
  }

  # Run the tunnel-state model simulation
  out <- ce_sim(
    pop = pop,
    probs = probs,
    inp = inp,
    c_unit = rc_unit,
    e_unit = re_unit,
    c_init = rc_init,
    e_init = re_init,
    distn = FALSE,
    interv_names = interv_names
  )
}

# Health Economic Analysis Outputs ---

# Calculate standard Cost-Effectiveness stats
he <- bcea(
  eff = -out$e_total / n_total,
  cost = -out$c_total / n_total,
  ref = 4,                # Reference strategy is no_icd (4th strategy)
  interventions = c("Cox 5-year SCD risk >4%",
                    "Cox 5-year SCD risk >6%",
                    "Observed Strategy",
                    "No ICD Strategy")
)

save(he, device, model_type, out, n_total, file = here::here("data/bcea_data.RData"))
