# =============================================================================
# CE Markov Model Analysis Script
# Description: Compares ICD vs Non-ICD strategies using BUGS transition probs
#              Stratified by 5-year risk prediction (Cox model)
# Author: N Green, UCL (Updated for robustness)
# =============================================================================

library(dplyr)
library(purrr)
library(BCEA)
library(reshape2)
library(readr)
library(here)
library(ggplot2)

# Configuration ---

TIME_HORIZON  <- 30    # Max time (years)
DISCOUNT_RATE <- 0.035
SCENARIO_FILE <- here::here("inst/extdata/scenario_input_values.csv")
BUGS_FILE     <- here::here("data/jags_output_all_strat.RData")
TRANS_FILE    <- here::here("data/transdat.RData")

# Model Structure
STATE_NAMES <- c("hcm", "scd/shock", "death")
N_STATES    <- length(STATE_NAMES)
# Note: The model uses 2 * N_STATES internally (ICD arm + No-ICD arm)

#=======
# data
#=======

load(BUGS_FILE)

# Extract lambda array directly: [Sims, Interv, From, To]
lambda_samples <- mm1$BUGSoutput$sims.list$lambda
n_sim <- dim(lambda_samples)[1]

# Load Transition Data (for names)
load(TRANS_FILE)
interv_names <- names(transdat)
n_intervs <- length(interv_names)

# individual patient data
# original extract
data("ipd_risk")

# model states:
# 1. HCM with ICD
# 2. shock
# 3. all-cause mortality
# 4. HCM (no ICD)
# 5. SCD
# 6. all-cause mortality

# structure:
# (1.icd  2.shock 3.mortality)
# (4.~icd 5.scd   6.mortality)


#=========
# set-up
#=========

# start state populations per intervention group
init_counts <-
  inits_by_interv(
    ipd_risk,
    group_names = interv_names) |>
  mutate(state = "hcm")

# Create full grid of starting states (filling zeros where needed)
init_full <- expand.grid(
  risk = interv_names,
  status = c("ICD", "low_risk"),
  state = STATE_NAMES) |>
  arrange(risk, status) |>
  left_join(init_counts, by = c("risk", "status", "state")) |>
  tidyr::replace_na(list(n = 0))

# Calculate total sample size dynamically (Robustness fix: no hardcoded 3672)
n_sample <- sum(init_full$n)

# Format population array for simulation: [State, Time, Sim]
# We split by risk group to create the list input for ce_sim
pop_list <- list()
pop_split <- split(init_full$n, init_full$risk)

for (nm in interv_names) {
  # Initialize array with zeros
  pop_list[[nm]] <- array(0, dim = c(2 * N_STATES, TIME_HORIZON, n_sim))
  # Set time=1 (start) population
  # Replicating start pop across sims (deterministic start)
  pop_list[[nm]][, 1, ] <- pop_split[[nm]]
}

# 4. Process Probabilities (Vectorized)
# =============================================================================
# Goal: Create list of arrays [From, To, Time, Sim] for each intervention

probs_list <- list()

# Vectorized adjustment for Shock -> Death mortality
# We assume Shock (State 2) carries same mortality risk as HCM-ICD (State 1)
# lambda_samples indices: [Sim, Interv, From, To]

# 1. Extract Mortality Probs (State 1 -> 3) for ALL sims and ALL interventions
p_death_bg <- lambda_samples[,, 1, 3]

# 2. Apply to Shock State (State 2)
# State 2->3 (Death) gets background mortality
lambda_samples[,, 2, 3] <- p_death_bg
# State 2->1 (Back to HCM) gets the remainder
lambda_samples[,, 2, 1] <- 1 - p_death_bg

# 3. Reshape into BCEA format
for (j in seq_along(interv_names)) {
  nm <- interv_names[j]

  # Create container: [From, To, Time(1), Sim]
  # Note: 2*N_STATES because we model ICD(1-3) and NoICD(4-6) simultaneously
  probs_array <- array(NA, dim = c(2 * N_STATES, 2 * N_STATES, 1, n_sim))

  # Extract slice for this intervention: [Sim, From, To]
  interv_slice <- lambda_samples[, j, , ]

  # Transpose to [From, To, Sim] and plug into array
  probs_array[,, 1, ] <- aperm(interv_slice, c(2, 3, 1))

  probs_list[[nm]] <- probs_array
}

# =============================================================================
# 5. Run Simulations from CSV Input
# =============================================================================

# Load Scenario Inputs
inp_tab <- read_csv(SCENARIO_FILE, col_names = TRUE, show_col_types = FALSE)

# Ensure Scenario_Label exists. If missing, create default labels.
if(!"Scenario_Label" %in% colnames(inp_tab)) {
  inp_tab$Scenario_Label <- paste0("Scen", inp_tab$scenario)
}

res_list <- list()

cat("Starting Simulations for", nrow(inp_tab), "scenarios...\n")

# Loop through sensitivity analysis scenarios
for (i in seq_len(nrow(inp_tab))) {

  if (i %% 5 == 0) cat("Processing scenario", i, "...\n")

  inp <- inp_tab[i, ]
  inp$t_repl <- 1 # no. tunnel states for replacement

  # run simulation for this specific ICD Version/Scenario
  res_list[[inp$Scenario_Label]] <- ce_sim(
    pop = pop_list,
    probs = probs_list,
    inp = inp,
    c_unit = rc_unit,
    e_unit = re_unit,
    c_init = rc_init,
    e_init = re_init,
    distn = FALSE,   # costs fixed means, Probabilities are PSA samples
    interv_names = interv_names
  )
}

# Save raw results
save(res_list, file = here::here("data/res_dsa_all_strat.RData"))


# =============================================================================
# 6. GRAND ANALYSIS: Automated Comparison of ICD Versions
# =============================================================================

# 1. Identify the Comparison Candidates
# ------------------------------------------------------------
# We identify the unique Labels (e.g., "vticd", "sicd")
# We assume the FIRST row for each label is the "Baseline" we want to compare.
unique_labels <- unique(inp_tab$Scenario_Label)
rows_to_compare <- sapply(unique_labels, function(lbl) which(inp_tab$Scenario_Label == lbl)[1])

cat("Comparing Baselines for:", paste(unique_labels, collapse=", "), "\n")

# 2. Initialize Data Containers with "No ICD"
# ------------------------------------------------------------
# We take No ICD from the first scenario (Row 1)
base_row <- 1
idx_noicd_strat <- which(interv_names == "no_icd")

grand_e <- res_list[[base_row]]$e_total[, idx_noicd_strat, drop=FALSE]
grand_c <- res_list[[base_row]]$c_total[, idx_noicd_strat, drop=FALSE]
grand_labels <- c("No ICD")

# 3. Define Strategies to extract for each version
# ------------------------------------------------------------
strats_to_compare <- c("risk_over_4", "risk_over_6")

# 4. Loop through unique ICD Versions and stack them
# ------------------------------------------------------------
for (i in seq_along(unique_labels)) {

  label <- unique_labels[i]
  row_idx <- rows_to_compare[i]

  for (strat in strats_to_compare) {
    col_idx <- which(interv_names == strat)

    # Extract from the specific scenario row
    e_temp <- res_list[[row_idx]]$e_total[, col_idx, drop=FALSE]
    c_temp <- res_list[[row_idx]]$c_total[, col_idx, drop=FALSE]

    # Append
    grand_e <- cbind(grand_e, e_temp)
    grand_c <- cbind(grand_c, c_temp)

    # Label: e.g., "vticd (risk_over_6)"
    grand_labels <- c(grand_labels, paste0(label, " (", strat, ")"))
  }
}

# 5. Run BCEA on the Grand Frontier
# ------------------------------------------------------------
he_grand <- bcea(
  eff = -grand_e / n_sample,  # Benefits (flip sign if your function returns losses)
  cost = -grand_c / n_sample, # Costs (flip sign if your function returns negative costs)
  ref = 1,                    # Reference is "No ICD"
  interventions = grand_labels
)

# 6. Visualization
# ------------------------------------------------------------

# A. The Cost-Effectiveness Plane
ceplane.plot(he_grand, graph = "ggplot2", title = "Grand Frontier: All ICD Versions") +
  theme_minimal() +
  scale_color_viridis_d()

ggsave(here::here("images/grand_frontier_ceplane.png"), dpi=300)

# B. The Cost-Effectiveness Acceptability Curve (CEAC)
gg_grand_ceac <- ceac.plot(he_grand, graph = "ggplot2", pos=c(0,1), currency="£") +
  theme_minimal() +
  labs(title = "Decision Analysis: Choosing ICD Version & Threshold") +
  geom_vline(xintercept = c(20000, 30000), linetype = "dashed")

ggsave(here::here("images/grand_frontier_ceac.png"), gg_grand_ceac, dpi=300)

# C. Net Monetary Benefit (NMB) Comparison
wtp_threshold <- 30000
nmb_grand <- (he_grand$e * wtp_threshold) - he_grand$c

# Convert to data frame for ggplot
nmb_df <- melt(nmb_grand, varnames = c("Sim", "Strategy"), value.name = "NMB")

ggplot(nmb_df, aes(x = Strategy, y = NMB, fill = Strategy)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  labs(title = paste0("Net Monetary Benefit at £", wtp_threshold), y = "NMB (£)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave(here::here("images/nmb_comparison_boxplot.png"), dpi=300)

