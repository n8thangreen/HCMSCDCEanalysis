# Script to run combined TV-ICD and S-ICD cost-effectiveness analysis (PSA)
# Compares all 6 strategies on the same Cost-Effectiveness Plane and CEAC plots
# Reference strategy is TV-ICD: Observed
# Uses parameter distributions from inst/extdata/input_distns.csv and 30-year tunnel model

library(dplyr)
library(purrr)
library(BCEA)
library(readr)
library(ggplot2)

devtools::load_all(here::here())

rpoint <- function(n, x) x

# 1. Load transition counts and patient risk profiles
load(here::here("data/transdat.RData"))
data("ipd_risk")

# 2. Load the direct BUGS MCMC output
load(here::here("data/jags_output_all_strat_direct.RData"))
lambda_samples <- mm1$BUGSoutput$sims.list$lambda
n_sim <- dim(lambda_samples)[1]
interv_names <- names(transdat)
n_total <- nrow(ipd_risk)
J <- 30
n_tunnel <- 15
S <- 2 * n_tunnel + 4

# 3. Read input parameter distributions
inp_dist <- readr::read_csv(
  here::here("inst/extdata/input_distns.csv"),
  col_names = TRUE, show_col_types = FALSE
)

init_counts_all <- inits_by_interv(ipd_risk, group_names = interv_names) |>
  mutate(state = "hcm")

init_full <- expand.grid(
  risk = interv_names,
  status = c("ICD", "low_risk"),
  state = c("hcm", "scd/shock", "death")) |>
  arrange(risk, status) |>
  left_join(init_counts_all, by = c("risk", "status", "state")) |>
  tidyr::replace_na(list(n = 0))

pop0 <- split(init_full$n, init_full$risk)

pop <- list()
for (i in interv_names) {
  pop[[i]] <- array(0, dim = c(S, J, n_sim))
  pop[[i]][1, 1, ] <- pop0[[i]][1]                  # ICD starts in State 1
  pop[[i]][2 * n_tunnel + 2, 1, ] <- pop0[[i]][4]   # No-ICD starts in State 2T+2
}

alpha <- 0

# =============================================================================
# A. SIMULATE TV-ICD (Tunnel state model with battery replacement)
# =============================================================================
cat("Simulating TV-ICD (Tunnel Model)...\n")

inp_tv <- inp_dist %>%
  dplyr::filter(label == "tvicd")

# Weibull device-replacement probabilities by tunnel state (n_sim x T), sampled like the MCMC lambda
get_p_tv <- function(v) as.numeric(inp_tv$params[inp_tv$variable == v])
set.seed(42)
p_repl_tv <- sample_repl_probs(n_sim, get_p_tv("repl_median"), get_p_tv("repl_shape"),
                               median_se = get_p_tv("repl_median_se"), shape_se = get_p_tv("repl_shape_se"))

probs_tv <- list()
for (j in seq_along(interv_names)) {
  nm <- interv_names[j]
  probs_array <- array(NA, dim = c(S, S, 1, n_sim))
  for (i in 1:n_sim) {
    probs_array[,, 1, i] <- lambda_tunnel(lambda_samples[i, j, , ], alpha = alpha, p_repl = p_repl_tv[i, ])
  }
  probs_tv[[nm]] <- probs_array
}

set.seed(42)

out_tv <- ce_sim(
  pop = pop,
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
cat("Simulating S-ICD (Tunnel Model)...\n")

inp_sicd <- inp_dist %>%
  dplyr::filter(label == "sicd")

# Weibull device-replacement probabilities by tunnel state (n_sim x T), sampled like the MCMC lambda
get_p_sicd <- function(v) as.numeric(inp_sicd$params[inp_sicd$variable == v])
set.seed(42)
p_repl_sicd <- sample_repl_probs(n_sim, get_p_sicd("repl_median"), get_p_sicd("repl_shape"),
                                 median_se = get_p_sicd("repl_median_se"), shape_se = get_p_sicd("repl_shape_se"))

probs_sicd <- list()
for (j in seq_along(interv_names)) {
  nm <- interv_names[j]
  probs_array <- array(NA, dim = c(S, S, 1, n_sim))
  for (i in 1:n_sim) {
    probs_array[,, 1, i] <- lambda_tunnel(lambda_samples[i, j, , ], alpha = alpha, p_repl = p_repl_sicd[i, ])
  }
  probs_sicd[[nm]] <- probs_array
}

set.seed(42)

out_sicd <- ce_sim(
  pop = pop,
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

cols_tv <- c("risk_over_4", "risk_over_6", "icd")

eff_tv <- out_tv$e_total[, cols_tv]
cost_tv <- out_tv$c_total[, cols_tv]
colnames(eff_tv) <- c("TV-ICD: Cox > 4%", "TV-ICD: Cox > 6%", "TV-ICD: Observed")
colnames(cost_tv) <- c("TV-ICD: Cox > 4%", "TV-ICD: Cox > 6%", "TV-ICD: Observed")

eff_sicd <- out_sicd$e_total[, cols_tv]
cost_sicd <- out_sicd$c_total[, cols_tv]
colnames(eff_sicd) <- c("S-ICD: Cox > 4%", "S-ICD: Cox > 6%", "S-ICD: Observed")
colnames(cost_sicd) <- c("S-ICD: Cox > 4%", "S-ICD: Cox > 6%", "S-ICD: Observed")

grand_eff <- cbind(eff_tv, eff_sicd)
grand_cost <- cbind(cost_tv, cost_sicd)

# Reference is TV-ICD: Observed (3rd column in grand_eff)
he_grand <- bcea(
  eff = -grand_eff / n_total,
  cost = -grand_cost / n_total,
  ref = 3,
  interventions = colnames(grand_eff)
)

print(summary(he_grand))

abs_costs <- colMeans(grand_cost / n_total)
abs_qalys <- colMeans(grand_eff / n_total)
df_abs <- data.frame(
  Strategy = colnames(grand_eff),
  Cost = abs_costs,
  QALYs = abs_qalys
)
print(df_abs)

# Save Cost-Effectiveness Plane
png(filename = here::here("images/ce_plane_combined.png"), width = 680, height = 680, res = 120)
ceplane.plot(he_grand, graph = "base", title = "Cost-Effectiveness Plane (TV-ICD & S-ICD vs TV-ICD Observed)")
dev.off()

# Generate pairwise comparisons against TV-ICD: Observed
ceac_list <- list()
k_vals <- he_grand$k
strategies_to_compare <- colnames(grand_eff)[colnames(grand_eff) != "TV-ICD: Observed"]

for (strat in strategies_to_compare) {
  he_temp <- bcea(
    eff = -grand_eff[, c(strat, "TV-ICD: Observed")] / n_total,
    cost = -grand_cost[, c(strat, "TV-ICD: Observed")] / n_total,
    ref = 2,
    interventions = c(strat, "TV-ICD: Observed")
  )
  ceac_list[[strat]] <- he_temp$ceac[, 1]
}

df_ceac_combined <- data.frame()
for (strat in strategies_to_compare) {
  df_ceac_combined <- rbind(df_ceac_combined, data.frame(
    k = k_vals,
    Probability = ceac_list[[strat]],
    Strategy = strat
  ))
}

gg_ceac_combined <- ggplot(df_ceac_combined, aes(x = k, y = Probability, color = Strategy, linetype = Strategy)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey55") +
  geom_vline(xintercept = 20000, linetype = "dashed", color = "grey55") +
  geom_vline(xintercept = 30000, linetype = "dashed", color = "grey55") +
  labs(
    title = "Cost-Effectiveness Acceptability Curve (30-Year Tunnel Model)",
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


# =============================================================================
# D. SAVE PAIRWISE PSA SUMMARY vs TV-ICD: Observed (used in the write-up)
# =============================================================================
k_wtp <- 25000
ref_nm <- "TV-ICD: Observed"
psa_summary <- do.call(rbind, lapply(colnames(grand_eff), function(strat) {
  dC <- (grand_cost[, strat] - grand_cost[, ref_nm]) / n_total
  dE <- (grand_eff[, strat] - grand_eff[, ref_nm]) / n_total
  nb <- k_wtp * dE - dC
  data.frame(
    Strategy = strat,
    Mean_Cost = mean(grand_cost[, strat]) / n_total,
    Mean_QALY = mean(grand_eff[, strat]) / n_total,
    Delta_C = mean(dC),
    Delta_E = mean(dE),
    ICER = mean(dC) / mean(dE),
    EIB = mean(nb),
    Prob_CE = mean(nb > 0)
  )
}))
print(psa_summary)
readr::write_csv(psa_summary, here::here("data/combined_psa_summary.csv"))
