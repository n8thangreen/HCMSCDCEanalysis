# Unified Cost-Effectiveness Markov Simulation Runner
#
# Use this script to run the tunnel-state Markov model analysis for different
# devices (TVICD vs. SICD). Note: setting inp$t_repl = 1 and c_icd_repl = 0
# recovers the 6-state non-replacement special case.

library(dplyr)
library(purrr)
library(BCEA)
library(readr)
library(ggplot2)
devtools::load_all(here::here())

# Set your configs ---

device <- "TVICD"       # Options: "TVICD" (Transvenous) or "SICD" (Subcutaneous)
J <- 30                 # Time horizon / cycle limit (years)

# 1. Load transition counts and patient risk profiles
data("transdat")
data("ipd_risk")

# 2. Load the BUGS MCMC output
load(here::here("data/jags_output.RData"))
lambda0_local <- mm1$BUGSoutput$sims.list$lambda.0
lambda1_local <- mm1$BUGSoutput$sims.list$lambda.1
lambda2_local <- mm1$BUGSoutput$sims.list$lambda.2
lambda3_local <- mm1$BUGSoutput$sims.list$lambda.3

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

# Setup initial starting populations from individual patient data
init_counts <- get_init_counts(ipd_risk)

# Tunnel-State Model Setup ---
n_tunnel <- inp$t_repl
S <- 2 * n_tunnel + 4

# Create starting population arrays
pop <- create_pop_array(init_counts, J, n_sim, n_tunnel = n_tunnel, interv_names = interv_names)

# Create transition probability arrays
alpha <- 0
lambda_list <- list(obs = lambda0_local, risk6 = lambda1_local, risk4 = lambda2_local, no_icd = lambda3_local)
probs <- create_probs_array(lambda_list, n_tunnel = n_tunnel, alpha = alpha)

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

# Health Economic Analysis Outputs ---

he <- bcea(
  eff = -out$e_total / n_total,
  cost = -out$c_total / n_total,
  ref = 3,                # Reference strategy obs
  interventions = c("Cox 5-year SCD risk >4%",
                    "Cox 5-year SCD risk >6%",
                    "Observed Strategy",
                    "No ICD Strategy")
)

model_type <- "tunnel"
save(he, device, model_type, out, n_total, file = here::here("data/bcea_data.RData"))
