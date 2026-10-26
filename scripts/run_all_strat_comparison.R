# =============================================================================
# CE Markov Model Analysis Script (Tunnel-State Replacement Model)
# Description: Compares Baseline TV-ICD (10-yr battery) vs Baseline S-ICD (8-yr battery)
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
WTP_THRESHOLD <- 25000
SCENARIO_FILE <- here::here("inst/extdata/scenario_input_values.csv")
BUGS_FILE     <- here::here("data/jags_output_all_strat.RData")
TRANS_FILE    <- here::here("data/transdat.RData")

# =============================================================================
# 1. LOAD DATA & SETUP POPULATIONS
# =============================================================================
cat("Loading Data...\n")

load(BUGS_FILE)
lambda_samples <- mm1$BUGSoutput$sims.list$lambda
n_sim <- dim(lambda_samples)[1]

load(TRANS_FILE)
interv_names <- names(transdat)
n_intervs <- length(interv_names)

data("ipd_risk")

init_counts <- inits_by_interv(ipd_risk, group_names = interv_names) |>
  mutate(state = "hcm")

init_full <- expand.grid(
  risk = interv_names,
  status = c("ICD", "low_risk"),
  state = c("hcm", "scd/shock", "death")) |>
  arrange(risk, status) |>
  left_join(init_counts, by = c("risk", "status", "state")) |>
  tidyr::replace_na(list(n = 0))

n_sample <- nrow(ipd_risk)
pop0 <- split(init_full$n, init_full$risk)

# =============================================================================
# 2. RUN TUNNEL SIMULATIONS (Baseline Rows for TV-ICD and S-ICD)
# =============================================================================

inp_tab <- read_csv(SCENARIO_FILE, col_names = TRUE, show_col_types = FALSE)

labels_found <- unique(inp_tab$Scenario_Label)
rows_to_run <- sapply(labels_found, function(x) which(inp_tab$Scenario_Label == x)[1])

cat("Comparing Baseline Rows:", paste(rows_to_run, collapse = ", "), "\n")

n_tunnel <- 15
S <- 2 * n_tunnel + 4
res_list <- list()

for (i in rows_to_run) {
  inp <- inp_tab[i, ]
  p_repl <- sample_repl_probs(n_sim, inp$repl_median, inp$repl_shape)  # deterministic Weibull (no SEs)
  cat("  -> Simulating Input Scenario Row", i, "(", inp$Scenario_Label, ")...\n")

  pop_list <- list()
  probs_list <- list()
  for (j in seq_along(interv_names)) {
    nm <- interv_names[j]
    pop_list[[nm]] <- array(0, dim = c(S, TIME_HORIZON, n_sim))
    pop_list[[nm]][1, 1, ] <- pop0[[nm]][1]
    pop_list[[nm]][2 * n_tunnel + 2, 1, ] <- pop0[[nm]][4]

    probs_array <- array(NA, dim = c(S, S, 1, n_sim))
    for (s_idx in 1:n_sim) {
      probs_array[,, 1, s_idx] <- lambda_tunnel(lambda_samples[s_idx, j, , ], alpha = 0, p_repl = p_repl[s_idx, ])
    }
    probs_list[[nm]] <- probs_array
  }

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
# 3. COMPARATIVE ANALYSIS: S-ICD vs TV-ICD across 8 Strategies
# =============================================================================

idx_tvicd_base <- rows_to_run[1]
idx_sicd_base  <- rows_to_run[2]

ce_data <- data.frame()

for (k in 1:n_intervs) {
  strat_name <- interv_names[k]

  eff_v <- res_list[[idx_tvicd_base]]$e_total[, k] / n_sample
  cost_v <- res_list[[idx_tvicd_base]]$c_total[, k] / n_sample

  eff_s <- res_list[[idx_sicd_base]]$e_total[, k] / n_sample
  cost_s <- res_list[[idx_sicd_base]]$c_total[, k] / n_sample

  # Incremental (S-ICD - TV-ICD)
  delta_e <- eff_s - eff_v
  delta_c <- cost_s - cost_v
  inmb <- (delta_e * WTP_THRESHOLD) - delta_c

  clean_label <- gsub("risk_over_", "Risk > ", strat_name)
  clean_label <- gsub("_", " ", clean_label)

  ce_data <- rbind(ce_data, data.frame(
    Strategy_ID = k,
    Description = clean_label,
    TVICD_Cost = mean(cost_v),
    TVICD_QALY = mean(eff_v),
    SICD_Cost = mean(cost_s),
    SICD_QALY = mean(eff_s),
    Delta_E_Mean = mean(delta_e),
    Delta_C_Mean = mean(delta_c),
    INMB_Mean = mean(inmb),
    INMB_Lower = quantile(inmb, 0.025),
    INMB_Upper = quantile(inmb, 0.975)
  ))
}

print(ce_data, row.names = FALSE)

save(ce_data, WTP_THRESHOLD, file = here::here("data/ce_analysis_all_strat_s-icd_results.RData"))
