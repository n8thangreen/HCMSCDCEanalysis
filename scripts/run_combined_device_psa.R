# Script to run combined TV-ICD and S-ICD cost-effectiveness analysis
# Compares all 6 strategies on the same Cost-Effectiveness Plane and CEAC plots
# Reference strategy is TV-ICD: Observed

library(dplyr)
library(purrr)
library(BCEA)
library(readr)
library(ggplot2)

devtools::load_all(here::here())

# 1. Load transition counts and patient risk profiles
load(here::here("data/transdat.RData"))
data("ipd_risk")

# 2. Load the direct BUGS MCMC output
load(here::here("data/jags_output_all_strat_direct.RData"))
lambda_samples <- mm1$BUGSoutput$sims.list$lambda
n_sim <- dim(lambda_samples)[1]
interv_names <- names(transdat)
n_total <- nrow(ipd_risk)
J <- 12

# Copy background mortality (ICD -> Death) to Shock state (Shock -> Death)
p_death_bg <- lambda_samples[,, 1, 3]
lambda_samples[,, 2, 3] <- p_death_bg
lambda_samples[,, 2, 1] <- 1 - p_death_bg

# 3. Read input parameter distributions
inp_dist <- readr::read_csv(
  here::here("inst/extdata/input_distns.csv"),
  col_names = TRUE, show_col_types = FALSE
)

# =============================================================================
# A. SIMULATE TV-ICD (Standard 6-state model)
# =============================================================================
cat("Simulating TV-ICD...\n")
S_tv <- 6
probs_tv <- list()
for (j in seq_along(interv_names)) {
  nm <- interv_names[j]
  probs_array <- array(NA, dim = c(2 * 3, 2 * 3, 1, n_sim))
  interv_slice <- lambda_samples[, j, , ]
  probs_array[,, 1, ] <- aperm(interv_slice, c(2, 3, 1))
  probs_tv[[nm]] <- probs_array
}

init_counts_tv <- inits_by_interv(ipd_risk, group_names = interv_names) |>
  mutate(state = "hcm")

init_full_tv <- expand.grid(
  risk = interv_names,
  status = c("ICD", "low_risk"),
  state = c("hcm", "scd/shock", "death")) |>
  arrange(risk, status) |>
  left_join(init_counts_tv, by = c("risk", "status", "state")) |>
  tidyr::replace_na(list(n = 0))

pop0_tv <- split(init_full_tv$n, init_full_tv$risk)
pop_tv <- list()
for (i in interv_names) {
  pop_tv[[i]] <- array(NA, dim = c(S_tv, J, n_sim))
  pop_tv[[i]][, 1, ] <- pop0_tv[[i]] # pop0_tv[[i]] is the 6-element vector: c(ICD, 0, 0, no_ICD, 0, 0)
}

# TV-ICD Parameter setup
inp_tv <- inp_dist %>%
  dplyr::filter(label == "vticd")

inp_tv <- inp_tv %>%
  mutate(params = case_when(
    variable == "c_appt" ~ "290",       # Point at 290
    variable == "p_compl" ~ "0.043",    # Point at 0.043
    variable == "t_repl" ~ "1",         # Set t_repl to 1 (no tunnel states)
    variable == "c_icd_repl" ~ "0",     # Set replacement cost to 0 (no replacements in standard model)
    variable == "q_hcm" ~ "0.88",       # Fix health parameters to match Table 2 QALYs
    variable == "u_icd" ~ "0.9",
    variable == "u_shock" ~ "0.875",
    variable == "u_implant" ~ "-0.048",
    TRUE ~ params
  )) %>%
  mutate(distn = case_when(
    variable %in% c("q_hcm", "u_icd", "u_shock", "u_implant") ~ "point",
    TRUE ~ distn
  ))

rpoint <- function(n, x) x
set.seed(42)

out_tv <- ce_sim(
  pop = pop_tv,
  probs = probs_tv,
  inp = inp_tv,
  c_unit = rc_unit,
  e_unit = re_unit,
  c_init = rc_init,
  e_init = re_init,
  distn = TRUE, # PSA
  interv_names = interv_names
)

# =============================================================================
# B. SIMULATE S-ICD (Tunnel state model with battery replacement)
# =============================================================================
cat("Simulating S-ICD...\n")
inp_sicd <- inp_dist %>%
  dplyr::filter(label == "sicd")

# Enforce S-ICD parameters
inp_sicd <- inp_sicd %>%
  mutate(params = case_when(
    variable == "c_appt" ~ "290",
    variable == "p_compl" ~ "0.047",
    variable == "t_repl" ~ "6",         # S-ICD battery lifetime is 6 years
    variable == "c_icd_repl" ~ "12000",
    variable == "c_icd" ~ "13788",      # S-ICD device cost
    variable == "q_hcm" ~ "0.88",
    variable == "u_icd" ~ "0.9",
    variable == "u_shock" ~ "0.875",
    variable == "u_implant" ~ "-0.048",
    TRUE ~ params
  )) %>%
  mutate(distn = case_when(
    variable %in% c("q_hcm", "u_icd", "u_shock", "u_implant", "c_icd", "c_icd_repl") ~ "point",
    TRUE ~ distn
  ))

n_tunnel <- 6
S_sicd <- 2 * n_tunnel + 4

probs_sicd <- list()
alpha <- 0
for (j in seq_along(interv_names)) {
  nm <- interv_names[j]
  probs_array <- array(NA, dim = c(S_sicd, S_sicd, 1, n_sim))
  for (i in 1:n_sim) {
    probs_array[,, 1, i] <- lambda_tunnel(lambda_samples[i, j, , ], n_tunnel = n_tunnel, alpha)
  }
  probs_sicd[[nm]] <- probs_array
}

pop_sicd <- list()
for (i in interv_names) {
  pop_sicd[[i]] <- array(NA, dim = c(S_sicd, J, n_sim))
  pop_sicd[[i]][, 1, ] <- 0
  pop_sicd[[i]][1, 1, ] <- pop0_tv[[i]][1] # ICD starts in State 1
  pop_sicd[[i]][2 * n_tunnel + 2, 1, ] <- pop0_tv[[i]][4] # No-ICD starts in State 2T+2 (4th element of pop0_tv)
}

out_sicd <- ce_sim(
  pop = pop_sicd,
  probs = probs_sicd,
  inp = inp_sicd,
  c_unit = rc_unit,
  e_unit = re_unit,
  c_init = rc_init,
  e_init = re_init,
  distn = TRUE, # PSA
  interv_names = interv_names
)

# =============================================================================
# C. COMBINE & ANALYZE
# =============================================================================
cat("Combining results...\n")
cols_tv <- c("risk_over_4", "risk_over_6", "icd")

# Select the 3 strategies for TV-ICD
eff_tv <- out_tv$e_total[, cols_tv]
cost_tv <- out_tv$c_total[, cols_tv]
colnames(eff_tv) <- c("TV-ICD: Cox > 4%", "TV-ICD: Cox > 6%", "TV-ICD: Observed")
colnames(cost_tv) <- c("TV-ICD: Cox > 4%", "TV-ICD: Cox > 6%", "TV-ICD: Observed")

# Select the 3 strategies for S-ICD
eff_sicd <- out_sicd$e_total[, cols_tv]
cost_sicd <- out_sicd$c_total[, cols_tv]
colnames(eff_sicd) <- c("S-ICD: Cox > 4%", "S-ICD: Cox > 6%", "S-ICD: Observed")
colnames(cost_sicd) <- c("S-ICD: Cox > 4%", "S-ICD: Cox > 6%", "S-ICD: Observed")

# Combine all 6 columns
grand_eff <- cbind(eff_tv, eff_sicd)
grand_cost <- cbind(cost_tv, cost_sicd)

# Reference is TV-ICD: Observed (which is 3rd column in grand_eff)
he_grand <- bcea(
  eff = -grand_eff / n_total,
  cost = -grand_cost / n_total,
  ref = 3, # TV-ICD: Observed is reference
  interventions = colnames(grand_eff)
)

# 1. Print summaries
cat("\n========================================================\n")
cat("=== COMBINED TV-ICD & S-ICD PSA COST-EFFECTIVENESS ===\n")
cat("========================================================\n")
print(summary(he_grand))

abs_costs <- colMeans(grand_cost / n_total)
abs_qalys <- colMeans(grand_eff / n_total)
df_abs <- data.frame(
  Strategy = colnames(grand_eff),
  Cost = abs_costs,
  QALYs = abs_qalys
)
print(df_abs)

# 2. Save Cost-Effectiveness Plane (using base R for plotting all 5 comparators relative to reference)
cat("\nGenerating Cost-effectiveness Plane plot (base R)...\n")
png(filename = here::here("images/ce_plane_combined.png"), width = 680, height = 680, res = 120)
ceplane.plot(he_grand, graph = "base", title = "Cost-Effectiveness Plane (TV-ICD & S-ICD vs TV-ICD Observed)")
dev.off()

# 3. Generate pairwise comparisons to plot pairwise CEAC curves against TV-ICD: Observed
cat("Generating CEAC plots...\n")
ceac_list <- list()
k_vals <- he_grand$k

# Compare each strategy to TV-ICD: Observed
strategies_to_compare <- colnames(grand_eff)[colnames(grand_eff) != "TV-ICD: Observed"]

for (strat in strategies_to_compare) {
  he_temp <- bcea(
    eff = -grand_eff[, c(strat, "TV-ICD: Observed")] / n_total,
    cost = -grand_cost[, c(strat, "TV-ICD: Observed")] / n_total,
    ref = 2,
    interventions = c(strat, "TV-ICD: Observed")
  )
  ceac_list[[strat]] <- he_temp$ceac[, 1] # Probability that comparator is cost-effective vs reference
}

# Combine CEAC probabilities into data frame
df_ceac_combined <- data.frame()
for (strat in strategies_to_compare) {
  df_ceac_combined <- rbind(df_ceac_combined, data.frame(
    k = k_vals,
    Probability = ceac_list[[strat]],
    Strategy = strat
  ))
}

# Plot combined CEAC curves
gg_ceac_combined <- ggplot(df_ceac_combined, aes(x = k, y = Probability, color = Strategy, linetype = Strategy)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey55") +
  geom_vline(xintercept = 20000, linetype = "dashed", color = "grey55") +
  geom_vline(xintercept = 30000, linetype = "dashed", color = "grey55") +
  labs(
    title = "Cost-Effectiveness Acceptability Curve",
    subtitle = "All strategies vs. TV-ICD: Observed",
    x = "Willingness-to-pay threshold k (£)",
    y = "Probability of cost-effectiveness"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(face = "italic", hjust = 0.5),
    legend.position = "bottom"
  ) +
  scale_color_brewer(palette = "Set1")

ggsave(
  filename = here::here("images/ceac_combined.png"),
  plot = gg_ceac_combined,
  width = 8,
  height = 7,
  dpi = 320
)

# Copy plots to artifact directory
dir.create(here::here("C:/Users/Nathan/.gemini/antigravity-cli/brain/c54c98ba-5eee-4e9b-a561-2dcb7a6deda4"), showWarnings = FALSE, recursive = TRUE)
file.copy(from = here::here("images/ce_plane_combined.png"), 
          to = "C:/Users/Nathan/.gemini/antigravity-cli/brain/c54c98ba-5eee-4e9b-a561-2dcb7a6deda4/ce_plane_combined.png", overwrite = TRUE)
file.copy(from = here::here("images/ceac_combined.png"), 
          to = "C:/Users/Nathan/.gemini/antigravity-cli/brain/c54c98ba-5eee-4e9b-a561-2dcb7a6deda4/ceac_combined.png", overwrite = TRUE)

cat("\nCombined analysis completed successfully! Plots saved to images/ and artifacts directory.\n")
