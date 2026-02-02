# =============================================================================
# CE Markov Model Analysis: Tunnel State Version
# Description: Compares ICD vs Non-ICD with Tunnel States for Device Replacement
#              (Tracks time-since-implant to apply replacement costs)
# Author: N Green, UCL (Updated for new Data Structure)
# =============================================================================

library(dplyr)
library(purrr)
library(BCEA)
library(reshape2)
library(readr)
library(here)
library(ggplot2)

# =============================================================================
# 1. CONFIGURATION & HELPER FUNCTIONS
# =============================================================================

TIME_HORIZON  <- 30
DISCOUNT_RATE <- 0.035
WTP_THRESHOLD <- 30000

# File Paths
SCENARIO_FILE <- here::here("scenario_input_values.csv")
BUGS_FILE     <- here::here("data/jags_output_all_strat.RData")
TRANS_FILE    <- here::here("data/transdat.RData")

# -----------------------------------------------------------------------------
# Helper: Create Tunnel Matrix
# Expands the 6x6 matrix into (2*T + 4)x(2*T + 4)
# -----------------------------------------------------------------------------
expand_to_tunnel <- function(mat_6x6, n_tunnel) {

  # New dimensions
  # ICD Arm: T states (Tunnel) + 1 Shock + 1 Death
  # NoICD Arm: T states (Tunnel) + 1 SCD + 1 Death
  # Total = 2*T + 4
  S_new <- 2 * n_tunnel + 4
  M <- matrix(0, nrow = S_new, ncol = S_new)

  # --- INDICES ---
  # ICD Arm
  idx_icd_start <- 1
  idx_shock     <- n_tunnel + 1
  idx_death_icd <- n_tunnel + 2

  # No-ICD Arm
  idx_noicd_start <- n_tunnel + 3
  idx_scd         <- 2 * n_tunnel + 3
  idx_death_noicd <- 2 * n_tunnel + 4

  # --- FILL ICD ARM ---
  # Transition probs from 6x6 matrix (Row 1)
  p_stay  <- mat_6x6[1, 1]
  p_shock <- mat_6x6[1, 2]
  p_death <- mat_6x6[1, 3]

  for (k in 1:n_tunnel) {
    curr <- idx_icd_start + k - 1

    # 1. Move to Next Tunnel State (Aging)
    # If at last tunnel state, loop back to start (Replacement Cycle)
    if (k < n_tunnel) {
      M[curr, curr + 1] <- p_stay
    } else {
      M[curr, idx_icd_start] <- p_stay # Replacement happens here
    }

    # 2. Risk Events (Shock/Death)
    M[curr, idx_shock]     <- p_shock
    M[curr, idx_death_icd] <- p_death
  }

  # Shock State (Fixed return to ICD Start)
  # Original assumption: Shock -> ICD (Next Cycle)
  # We send them to the first tunnel state (fresh start or continue?)
  # Usually: Return to same tunnel or specific state. Here mapped to Start.
  M[idx_shock, idx_icd_start] <- mat_6x6[2, 1]
  M[idx_shock, idx_death_icd] <- mat_6x6[2, 3] # Background mortality

  # Death State (Absorbent)
  M[idx_death_icd, idx_death_icd] <- 1


  # --- FILL NO-ICD ARM ---
  # Transition probs from 6x6 matrix (Row 4)
  p_stay_no  <- mat_6x6[4, 4]
  p_scd_no   <- mat_6x6[4, 5]
  p_death_no <- mat_6x6[4, 6]

  for (k in 1:n_tunnel) {
    curr <- idx_noicd_start + k - 1

    # 1. Move to Next (Symmetry with ICD arm, though no replacement cost)
    if (k < n_tunnel) {
      M[curr, curr + 1] <- p_stay_no
    } else {
      M[curr, idx_noicd_start] <- p_stay_no
    }

    # 2. Risk Events
    M[curr, idx_scd]         <- p_scd_no
    M[curr, idx_death_noicd] <- p_death_no
  }

  # SCD & Death (Absorbent)
  M[idx_scd, idx_scd] <- 1
  M[idx_death_noicd, idx_death_noicd] <- 1

  return(M)
}


# =============================================================================
# 2. LOAD DATA
# =============================================================================
cat("Loading Data...\n")

# BUGS Output
load(BUGS_FILE)
lambda_samples <- mm1$BUGSoutput$sims.list$lambda # [Sim, Interv, From, To]
n_sim <- dim(lambda_samples)[1]

# Strategy Names
load(TRANS_FILE)
interv_names <- names(transdat)
n_intervs <- length(interv_names)

# Patient Data
data("ipd_risk")

# Scenario Inputs (for n_tunnel)
inp_tab <- read_csv(SCENARIO_FILE, col_names = TRUE, show_col_types = FALSE)
if(!"Scenario_Label" %in% colnames(inp_tab)) inp_tab$Scenario_Label <- paste0("Scen_", 1:nrow(inp_tab))

# Get Tunnel Count (Time to Replacement) from first row
# (Assuming it's constant for the structure)
n_tunnel <- inp_tab$t_repl[1]
S_tunnel <- 2 * n_tunnel + 4
cat("Tunnel States configured. n_tunnel =", n_tunnel, "| Total States =", S_tunnel, "\n")


# =============================================================================
# 3. SETUP POPULATIONS (Mapped to Tunnel)
# =============================================================================
cat("Setting up Populations...\n")

# Calculate Initial Counts (Total)
init_counts <- inits_by_interv(ipd_risk, group_names = interv_names) |>
  mutate(state = "hcm")

init_full <- expand.grid(risk = interv_names, status = c("ICD", "low_risk"), state = STATE_NAMES) |>
  arrange(risk, status) |>
  left_join(init_counts, by = c("risk", "status", "state")) |>
  tidyr::replace_na(list(n = 0))

n_sample <- sum(init_full$n)

# Format into Tunnel Array: [State, Time, Sim]
pop_list <- list()
pop_split <- split(init_full$n, init_full$risk)

# Indicies for Start States in the Tunnel Matrix
idx_icd_start   <- 1
idx_noicd_start <- n_tunnel + 3

for (nm in interv_names) {
  # Create Empty Tunnel Population Array
  pop_arr <- array(0, dim = c(S_tunnel, TIME_HORIZON, n_sim))

  # Get counts for this strategy
  # Standard vector order: [ICD, Shock, Death, NoICD, SCD, Death]
  counts <- pop_split[[nm]]

  # Map ICD count -> First ICD Tunnel State
  pop_arr[idx_icd_start, 1, ] <- counts[1] # ICD Healthy

  # Map No-ICD count -> First No-ICD Tunnel State
  pop_arr[idx_noicd_start, 1, ] <- counts[4] # NoICD Healthy

  # Note: Usually we start with 0 in Shock/Death, so other mappings aren't needed
  # unless you have prevalent cases.

  pop_list[[nm]] <- pop_arr
}


# =============================================================================
# 4. PROCESS PROBABILITIES (Expanding to Tunnel)
# =============================================================================
cat("Processing Probabilities & Creating Tunnels...\n")

probs_list <- list()

# Pre-process background mortality for Shock
p_death_bg <- lambda_samples[,, 1, 3]
lambda_samples[,, 2, 3] <- p_death_bg
lambda_samples[,, 2, 1] <- 1 - p_death_bg

# Loop Strategies
for (j in seq_along(interv_names)) {
  nm <- interv_names[j]

  # Container: [From, To, Time, Sim] with NEW Dimensions
  probs_array <- array(NA, dim = c(S_tunnel, S_tunnel, 1, n_sim))

  # Loop Simulations
  for (i in 1:n_sim) {
    # Extract the 6x6 matrix for this Sim & Strat
    mat_6x6 <- lambda_samples[i, j, , ]

    # Expand to Tunnel Matrix
    mat_tunnel <- expand_to_tunnel(mat_6x6, n_tunnel)

    # Store
    probs_array[,, 1, i] <- mat_tunnel
  }

  probs_list[[nm]] <- probs_array
}


# =============================================================================
# 5. RUN SIMULATIONS
# =============================================================================
# We assume you want to run the 8 comparison scenarios as discussed previously
# (Rows 1-8 VTICD vs 19-26 SICD)

rows_vticd <- 1:8
rows_sicd  <- 19:26
rows_to_run <- c(rows_vticd, rows_sicd)

res_list <- list()

cat("Running Tunnel Simulations for", length(rows_to_run), "Scenarios...\n")

for (i in rows_to_run) {

  cat("  -> Simulating Row", i, "...\n")

  inp <- inp_tab[i, ]

  # Ensure t_repl in input matches config
  inp$t_repl <- n_tunnel

  res_list[[i]] <- ce_sim(
    pop = pop_list,
    probs = probs_list,
    inp = inp,
    c_unit = rc_unit,
    e_unit = re_unit,
    c_init = rc_init,
    e_init = re_init,
    distn = FALSE,
    interv_names = interv_names
  )
}


# =============================================================================
# 6. GRAND ANALYSIS (INMB Forest Plot)
# =============================================================================
cat("Generating Comparisons...\n")

ce_data <- data.frame()
thresholds <- c("risk_over_4", "risk_over_6")

# Descriptions for the 8 scenarios
scen_descriptions <- c("Baseline", "Low QoL", "High QoL", "High Utility Shock",
                       "No Disutility", "Low Battery", "Perfect Utility", "Low Cost")

for (k in 1:8) {
  idx_v <- rows_vticd[k]
  idx_s <- rows_sicd[k]
  desc  <- scen_descriptions[k]

  for (strat in thresholds) {
    col_idx <- which(interv_names == strat)

    # Extract Means
    # NOTE: Check if ce_sim output needs sign flip (-ve costs)
    eff_v <- mean(-res_list[[idx_v]]$e_total[, col_idx]) / n_sample
    cost_v <- mean(-res_list[[idx_v]]$c_total[, col_idx]) / n_sample

    eff_s <- mean(-res_list[[idx_s]]$e_total[, col_idx]) / n_sample
    cost_s <- mean(-res_list[[idx_s]]$c_total[, col_idx]) / n_sample

    # Calculate INMB (SICD - VTICD)
    delta_e <- eff_s - eff_v
    delta_c <- cost_s - cost_v
    inmb <- (delta_e * WTP_THRESHOLD) - delta_c

    ce_data <- rbind(ce_data, data.frame(
      Scenario = desc,
      Threshold = gsub("risk_over_", "> ", strat),
      INMB = inmb
    ))
  }
}

# Plot
gg_forest <- ggplot(ce_data, aes(x = INMB, y = reorder(Scenario, INMB), color = Threshold)) +
  geom_vline(xintercept = 0, color = "black") +
  geom_point(size = 4) +
  labs(title = "Incremental Net Benefit (Tunnel Model): SICD vs VTICD",
       x = "INMB (£)", y = "Scenario") +
  theme_minimal()

ggsave(here::here("images/tunnel_inmb_forest.png"), gg_forest, width = 8, height = 6)

cat("Tunnel Analysis Complete.\n")
