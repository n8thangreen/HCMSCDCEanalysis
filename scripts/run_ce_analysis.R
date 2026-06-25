# Unified Cost-Effectiveness Markov Simulation Runner
#
# Use this script to run the analysis for different devices (ICD vs. SICD)
# and model structures (Standard vs. Tunnel state models).

library(dplyr)
library(purrr)
library(BCEA)
library(readr)
library(ggplot2)
library(HCMSCDCEanalysis)

# Set your configs ---

device <- "SICD"       # Options: "ICD" (Transvenous) or "SICD" (Subcutaneous)
model_type <- "standard" # Options: "standard" (6-state model) or "tunnel" (with replacement cycles)
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
lbl <- if (tolower(device) == "sicd") "sicd" else "tvicd"
inp <- inp_tab %>%
  dplyr::filter(Scenario_Label == lbl) %>%
  dplyr::slice(1)

# Model Structuring (Standard vs. Tunnel) ---

if (model_type == "standard") {
  S <- 6  # 6-state model
  inp$t_repl <- 1 # Enforce no tunnel states for standard 6-state model

  # Starting populations
  init <- inits_by_interv(ipd_risk, group_names = c("icd", "risk_over_6", "risk_over_4"))
  n_init <- list(
    obs = c(filter(init, risk == "icd", status == "ICD")$n, 0, 0, filter(init, risk == "icd", status == "low_risk")$n, 0, 0),
    risk6 = c(filter(init, risk == "risk_over_6", status == "ICD")$n, 0, 0, filter(init, risk == "risk_over_6", status == "low_risk")$n, 0, 0),
    risk4 = c(filter(init, risk == "risk_over_4", status == "ICD")$n, 0, 0, filter(init, risk == "risk_over_4", status == "low_risk")$n, 0, 0),
    no_icd = c(0, 0, 0, n_total, 0, 0)
  )

  # Probability transition lists
  probs_empty <- array(NA, dim = c(S, S, 1, n_sim))
  probs <- list(obs = probs_empty, risk6 = probs_empty, risk4 = probs_empty, no_icd = probs_empty)

  for (i in 1:n_sim) {
    lambda.0[i, 2, 3] <- lambda.0[i, 1, 3]
    lambda.1[i, 2, 3] <- lambda.1[i, 1, 3]
    lambda.2[i, 2, 3] <- lambda.2[i, 1, 3]
    lambda.3[i, 2, 3] <- lambda.3[i, 1, 3]

    lambda.0[i, 2, 1] <- 1 - lambda.0[i, 2, 3]
    lambda.1[i, 2, 1] <- 1 - lambda.1[i, 2, 3]
    lambda.2[i, 2, 1] <- 1 - lambda.2[i, 2, 3]
    lambda.3[i, 2, 1] <- 1 - lambda.3[i, 2, 3]
  }

  probs$obs[,,1,] <- aperm(lambda.0, c(2,3,1))
  probs$risk6[,,1,] <- aperm(lambda.1, c(2,3,1))
  probs$risk4[,,1,] <- aperm(lambda.2, c(2,3,1))
  probs$no_icd[,,1,] <- aperm(lambda.3, c(2,3,1))

  # Populate starting cohort state populations
  pop <- list()
  for (i in interv_names) {
    pop[[i]] <- array(NA, dim = c(S, J, n_sim))
    pop[[i]][, 1, ] <- t(n_init[[i]])
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

  probs_empty <- array(NA, dim = c(S, S, 1, n_sim))
  probs <- list(
    obs = probs_empty,
    risk6 = probs_empty,
    risk4 = probs_empty,
    no_icd = probs_empty
  )

  alpha <- 0
  for (i in 1:n_sim) {
    probs$obs[,,1,i] <- lambda_tunnel(lambda.0[i,,], n_tunnel = n_tunnel, alpha)
    probs$risk6[,,1,i] <- lambda_tunnel(lambda.1[i,,], n_tunnel = n_tunnel, alpha)
    probs$risk4[,,1,i] <- lambda_tunnel(lambda.2[i,,], n_tunnel = n_tunnel, alpha)
    probs$no_icd[,,1,i] <- lambda_tunnel(lambda.3[i,,], n_tunnel = n_tunnel, alpha)
  }

  # Sample initial populations using risk score distributions
  # (Simulates patient cohorts entering tunnel replacement states)
  pop <- list()
  n_init_fixed <- c(
    risk4 = 2569,
    risk6 = 3134,
    obs = 3113,
    no_icd = 0
  )

  for (i in interv_names) {
    pop[[i]] <- array(NA, dim = c(S, J, n_sim))

    # Place initial patients in the first tunnel state (or healthy no-ICD state)
    pop[[i]][, 1, ] <- 0
    pop[[i]][1, 1, ] <- n_init_fixed[i]
    pop[[i]][n_tunnel + 1, 1, ] <- n_total - n_init_fixed[i]
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

# Health Economic Analysis (BCEA) and Outputs ---

# Calculate standard Cost-Effectiveness stats
he <- bcea(
  eff = -out$e_total / n_total,
  cost = -out$c_total / n_total,
  ref = 4, # Reference strategy is no_icd (4th strategy)
  interventions = c("Cox 5-year SCD risk >4%",
                    "Cox 5-year SCD risk >6%",
                    "Observed Strategy",
                    "No ICD Strategy")
)

# Plot and save output figures
dir.create(here::here("images"), showWarnings = FALSE)

# 1. Cost-effectiveness plane
ceplane.plot(he, currency = "£", graph = "ggplot", title = paste("CE Plane for", device))
ggsave(filename = paste0("images/ce_plane_", device, "_", model_type, ".png"), width = 7, height = 7, dpi = 320)

# 2. Cost-effectiveness acceptability curve (CEAC)
ceac.plot(he, title = paste("CEAC for", device), graph = "ggplot", currency = "£") +
  theme(text = element_text(size = 14)) +
  geom_hline(yintercept = 0.5, col = "grey", linetype = "dashed")

ggsave(filename = paste0("images/ceac_", device, "_", model_type, ".png"), width = 7, height = 7, dpi = 320)

# Print Summary Results
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("Device Type : ", device, "\n")
cat("Model Type  : ", model_type, "\n")
print(summary(he))
