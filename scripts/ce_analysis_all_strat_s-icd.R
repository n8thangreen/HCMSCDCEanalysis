# =============================================================================
# CE Markov Model Analysis Script
# Description: Compares Baseline VTICD vs Baseline SICD
#              across all 8 Risk Stratification (Transition) Scenarios
# Author: N Green, UCL
# =============================================================================

library(dplyr)
library(purrr)
library(BCEA)
library(reshape2)
library(readr)
library(here)
library(ggplot2)
library(ggrepel)
devtools::load_all(here::here())

# Configuration ---
TIME_HORIZON  <- 30
DISCOUNT_RATE <- 0.035
WTP_THRESHOLD <- 30000
SCENARIO_FILE <- here::here("inst/extdata/scenario_input_values.csv")
BUGS_FILE     <- here::here("data/jags_output_all_strat.RData")
TRANS_FILE    <- here::here("data/transdat.RData")

# Model Structure
STATE_NAMES <- c("hcm", "scd/shock", "death")
N_STATES    <- length(STATE_NAMES)

# =============================================================================
# 1. LOAD DATA & SETUP POPULATIONS
# =============================================================================
cat("Loading Data...\n")

# BUGS Output (Transition Probabilities)
load(BUGS_FILE)
lambda_samples <- mm1$BUGSoutput$sims.list$lambda
n_sim <- dim(lambda_samples)[1]

# Strategy Names (The "8 scenarios" related to transition probs)
load(TRANS_FILE)
interv_names <- names(transdat)
# e.g., "risk_over_4", "risk_over_6", "obs", "any_rf", etc.
n_intervs <- length(interv_names)

# Patient Data
data("ipd_risk")

# Setup Populations
init_counts <- inits_by_interv(ipd_risk, group_names = interv_names) |>
  mutate(state = "hcm")

init_full <- expand.grid(
  risk = interv_names,
  status = c("ICD", "low_risk"),
  state = STATE_NAMES) |>
  arrange(risk, status) |>
  left_join(init_counts, by = c("risk", "status", "state")) |>
  tidyr::replace_na(list(n = 0))

n_sample <- sum(init_full$n)

# Format population array
pop_list <- list()
pop_split <- split(init_full$n, init_full$risk)
for (nm in interv_names) {
  pop_list[[nm]] <- array(0, dim = c(2 * N_STATES, TIME_HORIZON, n_sim))
  pop_list[[nm]][, 1, ] <- pop_split[[nm]]
}

# =============================================================================
# 2. PROCESS PROBABILITIES
# =============================================================================
cat("Processing Probabilities...\n")
probs_list <- list()
p_death_bg <- lambda_samples[,, 1, 3]
lambda_samples[,, 2, 3] <- p_death_bg
lambda_samples[,, 2, 1] <- 1 - p_death_bg

for (j in seq_along(interv_names)) {
  nm <- interv_names[j]
  probs_array <- array(NA, dim = c(2 * N_STATES, 2 * N_STATES, 1, n_sim))
  interv_slice <- lambda_samples[, j, , ]
  probs_array[,, 1, ] <- aperm(interv_slice, c(2, 3, 1))
  probs_list[[nm]] <- probs_array
}

# =============================================================================
# 3. RUN SIMULATIONS (Baseline Rows Only)
# =============================================================================
# We only want the FIRST row of each Scenario_Label (vticd, sicd)
# This fixes the input parameters (Costs/Utilities)

inp_tab <- read_csv(SCENARIO_FILE, col_names = TRUE, show_col_types = FALSE)

# Robustness: Find the indices dynamically
labels_found <- unique(inp_tab$Scenario_Label)
rows_to_run <- sapply(labels_found, function(x) which(inp_tab$Scenario_Label == x)[1])

cat("Comparing Baseline Rows:", paste(rows_to_run, collapse=", "), "\n")
# Should be Row 1 (vticd) and Row 19 (sicd) based on your file

res_list <- list()

for (i in rows_to_run) {
  cat("  -> Simulating Input Scenario Row", i, "...\n")
  inp <- inp_tab[i, ]
  inp$t_repl <- 1

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
# 4. COMPARATIVE ANALYSIS: SICD vs VTICD across 8 Strategies
# =============================================================================

# Identify the Result Indices
# Assuming 'vticd' comes first in your CSV, and 'sicd' comes second
idx_vticd_base <- rows_to_run[1]
idx_sicd_base  <- rows_to_run[2]

ce_data <- data.frame()

# Loop through the 8 Transition Scenarios (Strategies)
for (k in 1:n_intervs) {

  strat_name <- interv_names[k]

  # Extract Results for this strategy
  # Note: -ve sign assumes ce_sim returns negative values for plotting conventions
  eff_v <- -res_list[[idx_vticd_base]]$e_total[, k] / n_sample
  cost_v <- -res_list[[idx_vticd_base]]$c_total[, k] / n_sample

  eff_s <- -res_list[[idx_sicd_base]]$e_total[, k] / n_sample
  cost_s <- -res_list[[idx_sicd_base]]$c_total[, k] / n_sample

  # Incremental (SICD - VTICD)
  delta_e <- eff_s - eff_v
  delta_c <- cost_s - cost_v
  inmb <- (delta_e * WTP_THRESHOLD) - delta_c

  # Clean up label (e.g. "risk_over_4" -> "> 4%")
  clean_label <- gsub("risk_over_", "Risk > ", strat_name)
  clean_label <- gsub("_", " ", clean_label)

  ce_data <- rbind(ce_data, data.frame(
    Strategy_ID = k,
    Description = clean_label,
    Delta_E_Mean = mean(delta_e),
    Delta_C_Mean = mean(delta_c),
    INMB_Mean = mean(inmb),
    # For error bars (optional)
    INMB_Lower = quantile(inmb, 0.025),
    INMB_Upper = quantile(inmb, 0.975)
  ))
}

# =============================================================================
# 5. VISUALIZATION
# =============================================================================

# A. INMB Forest Plot (The Decision Plot)
# -------------------
# This shows which risk strategy makes the SICD most cost-effective vs VTICD
gg_forest <- ggplot(ce_data, aes(x = INMB_Mean, y = reorder(Description, INMB_Mean))) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_errorbarh(aes(xmin = INMB_Lower, xmax = INMB_Upper), height = 0.2, color = "grey40") +
  geom_point(size = 4, color = "darkblue") +

  geom_text(aes(label = round(INMB_Mean, 0)), vjust = -1.2, size = 3.5) +

  labs(
    title = "Incremental Net Benefit: SICD vs VTICD",
    subtitle = "Comparing across 8 Risk Stratification Scenarios",
    x = paste0("INMB (£) at WTP £", WTP_THRESHOLD),
    y = "Risk Strategy",
    caption = "Positive values favor SICD. Error bars are 95% CrI."
  ) +
  theme_minimal()

gg_forest

ggsave(here::here("images/inmb_strategies_forest.png"), gg_forest, width = 8, height = 6, dpi = 300)


# B. Cost-Effectiveness Plane (Strategy Scatter)
# ----------------------------------------------
# Shows the Cost vs Benefit trade-off for each strategy
gg_plane <- ggplot(ce_data, aes(x = Delta_E_Mean, y = Delta_C_Mean, label = Description)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  geom_abline(intercept = 0, slope = WTP_THRESHOLD, linetype = "dashed", color = "grey60") +

  geom_point(size = 4, color = "darkred", alpha = 0.8) +
  geom_text_repel(size = 3.5) +

  labs(
    title = "Cost-Effectiveness Plane (SICD vs VTICD)",
    subtitle = "Mean outcomes for each Risk Strategy",
    x = "Incremental QALYs (SICD - VTICD)",
    y = "Incremental Cost (£) (SICD - VTICD)"
  ) +
  theme_minimal()

gg_plane

ggsave(here::here("images/ce_plane_strategies.png"), gg_plane, width = 8, height = 6, dpi = 300)

