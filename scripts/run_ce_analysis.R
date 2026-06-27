# Unified Cost-Effectiveness Markov Simulation Runner
#
# Use this script to run the analysis for different devices (ICD vs. SICD)
# and model structures (static vs. Tunnel state models).

library(dplyr)
library(purrr)
library(BCEA)
library(readr)
library(ggplot2)
devtools::load_all(here::here())

# Set your configs ---

device <- "TVICD"       # Options: "TVICD" (Transvenous) or "SICD" (Subcutaneous)
model_type <- "static" # Options: "static" (6-state model) or "tunnel" (with replacement cycles)
J <- 30                  # Time horizon / cycle limit (years)

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
lbl <- tolower(device)
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
init_counts <- get_init_counts(ipd_risk)

# Model Structuring (static vs. Tunnel) ---

if (model_type == "static") {
  S <- 6  # 6-state model
  inp$t_repl <- 1 # Enforce no tunnel states for static 6-state model

  # Create starting population arrays
  pop <- create_pop_array(init_counts, J, n_sim, model_type = "static", interv_names = interv_names)

  # Create transition probability arrays
  lambda_list <- list(obs = lambda.0, risk6 = lambda.1, risk4 = lambda.2, no_icd = lambda.3)
  probs <- create_probs_array(lambda_list, model_type = "static")

  # Run the static 6-state simulation
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

  # Create starting population arrays
  pop <- create_pop_array(init_counts, J, n_sim, model_type = "tunnel", n_tunnel = n_tunnel, interv_names = interv_names)

  # Create transition probability arrays
  alpha <- 0
  lambda_list <- list(obs = lambda.0, risk6 = lambda.1, risk4 = lambda.2, no_icd = lambda.3)
  probs <- create_probs_array(lambda_list, model_type = "tunnel", n_tunnel = n_tunnel, alpha = alpha)

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

# Calculate static Cost-Effectiveness stats
he <- bcea(
  eff = -out$e_total / n_total,
  cost = -out$c_total / n_total,
  ref = 3,                # Reference strategy obs
  interventions = c("Cox 5-year SCD risk >4%",
                    "Cox 5-year SCD risk >6%",
                    "Observed Strategy",
                    "No ICD Strategy")
)

save(he, device, model_type, out, n_total, file = here::here("data/bcea_data.RData"))
