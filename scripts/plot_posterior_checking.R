# Posterior Model Checking & Predictive Checks
# HCMSCDCEanalysis Package
#

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)
library(here)

# 1. Load data and MCMC output
load(here::here("data/transdat.RData"))
load(here::here("data/jags_output.RData"))

# 2. Extract raw counts and calculate sample sizes/proportions
data_obs <- transdat$icd
data_risk6 <- transdat$risk_over_6
data_risk4 <- transdat$risk_over_4
data_no_icd <- transdat$no_icd

n_S <- nrow(data_obs$ICD)
empty_mat <- matrix(0, nrow = n_S, ncol = n_S)

r.0 <- rbind(cbind(data_obs$ICD[, 1:n_S], empty_mat),
             cbind(empty_mat, data_obs$low_risk[, 1:n_S]))

r.1 <- rbind(cbind(data_risk6$ICD[, 1:n_S], empty_mat),
             cbind(empty_mat, data_risk6$low_risk[, 1:n_S]))

r.2 <- rbind(cbind(data_risk4$ICD[, 1:n_S], empty_mat),
             cbind(empty_mat, data_risk4$low_risk[, 1:n_S]))

r.3 <- rbind(cbind(data_no_icd$ICD[, 1:n_S], empty_mat),
             cbind(empty_mat, data_no_icd$low_risk[, 1:n_S]))

n.0 <- unname(rowSums(r.0))
n.1 <- unname(rowSums(r.1))
n.2 <- unname(rowSums(r.2))
n.3 <- unname(rowSums(r.3))

n_icd <- c(n.0[1], n.1[1], n.2[1], n.3[1])
n_noicd <- c(n.0[4], n.1[4], n.2[4], n.3[4])

r_icd_shock <- c(r.0[1, 2], r.1[1, 2], r.2[1, 2], r.3[1, 2])
r_icd_death <- c(r.0[1, 3], r.1[1, 3], r.2[1, 3], r.3[1, 3])
r_noicd_scd <- c(r.0[4, 5], r.1[4, 5], r.2[4, 5], r.3[4, 5])
r_noicd_death <- c(r.0[4, 6], r.1[4, 6], r.2[4, 6], r.3[4, 6])

strategies <- c("Observed", "Risk > 6%", "Risk > 4%", "No ICD")

raw_probs <- tibble(
  Strategy = rep(strategies, each = 4),
  State = rep(c("ICD", "ICD", "No ICD", "No ICD"), 4),
  Transition = rep(c("ICD -> Shock", "ICD -> Death", "No ICD -> SCD", "No ICD -> Death"), 4),
  Observed_Val = c(
    r_icd_shock[1] / n_icd[1], r_icd_death[1] / n_icd[1], r_noicd_scd[1] / n_noicd[1], r_noicd_death[1] / n_noicd[1],
    r_icd_shock[2] / n_icd[2], r_icd_death[2] / n_icd[2], r_noicd_scd[2] / n_noicd[2], r_noicd_death[2] / n_noicd[2],
    r_icd_shock[3] / n_icd[3], r_icd_death[3] / n_icd[3], r_noicd_scd[3] / n_noicd[3], r_noicd_death[3] / n_noicd[3],
    NaN, NaN, r_noicd_scd[4] / n_noicd[4], r_noicd_death[4] / n_noicd[4]
  )
) %>%
  filter(!is.nan(Observed_Val) & !is.na(Observed_Val))

raw_counts <- tibble(
  Strategy = rep(strategies, each = 4),
  Transition = rep(c("ICD -> Shock", "ICD -> Death", "No ICD -> SCD", "No ICD -> Death"), 4),
  Observed_Count = c(
    r_icd_shock[1], r_icd_death[1], r_noicd_scd[1], r_noicd_death[1],
    r_icd_shock[2], r_icd_death[2], r_noicd_scd[2], r_noicd_death[2],
    r_icd_shock[3], r_icd_death[3], r_noicd_scd[3], r_noicd_death[3],
    r_icd_shock[4], r_icd_death[4], r_noicd_scd[4], r_noicd_death[4]
  )
) %>%
  filter(!(Strategy == "No ICD" & Transition %in% c("ICD -> Shock", "ICD -> Death")))

# 3. Extract MCMC posterior draws
sims <- mm1$BUGSoutput$sims.list
n_sims <- mm1$BUGSoutput$n.sims

# Reconstruct posterior draws of transition probabilities
prob_draws_list <- list()
for (i in 1:4) {
  strat_name <- strategies[i]
  lambda_var <- paste0("lambda.", i - 1)
  lambda_draws <- sims[[lambda_var]] # shape [n_sims, 6, 6]

  df_icd_shock <- tibble(
    Strategy = strat_name,
    Transition = "ICD -> Shock",
    Value = lambda_draws[, 1, 2]
  )
  df_icd_death <- tibble(
    Strategy = strat_name,
    Transition = "ICD -> Death",
    Value = lambda_draws[, 1, 3]
  )
  df_noicd_scd <- tibble(
    Strategy = strat_name,
    Transition = "No ICD -> SCD",
    Value = lambda_draws[, 4, 5]
  )
  df_noicd_death <- tibble(
    Strategy = strat_name,
    Transition = "No ICD -> Death",
    Value = lambda_draws[, 4, 6]
  )

  prob_draws_list[[i]] <- bind_rows(df_icd_shock, df_icd_death, df_noicd_scd, df_noicd_death)
}

prob_draws <- bind_rows(prob_draws_list) %>%
  filter(!(Strategy == "No ICD" & Transition %in% c("ICD -> Shock", "ICD -> Death")))

# 4. Perform Posterior Predictive Check (PPC) simulation of transition counts
count_draws_list <- list()
for (i in 1:4) {
  strat_name <- strategies[i]

  # Simulate ICD state counts if there are patients in ICD state
  if (n_icd[i] > 0) {
    rates_icd_shock <- sims$rate_icd_shock[, i]
    rates_icd_death <- sims$rate_icd_death[, i]

    sim_icd_shock <- rbinom(n_sims, size = n_icd[i], prob = rates_icd_shock)
    sim_icd_death <- rbinom(n_sims, size = n_icd[i], prob = rates_icd_death)

    count_draws_list[[length(count_draws_list) + 1]] <- tibble(
      Strategy = strat_name,
      Transition = "ICD -> Shock",
      Value = sim_icd_shock
    )
    count_draws_list[[length(count_draws_list) + 1]] <- tibble(
      Strategy = strat_name,
      Transition = "ICD -> Death",
      Value = sim_icd_death
    )
  }

  # Simulate No ICD state counts (always > 0)
  rates_noicd_scd <- sims$rate_noicd_scd[, i]
  rates_noicd_death <- sims$rate_noicd_death[, i]

  sim_noicd_scd <- rbinom(n_sims, size = n_noicd[i], prob = rates_noicd_scd)
  sim_noicd_death <- rbinom(n_sims, size = n_noicd[i], prob = rates_noicd_death)

  count_draws_list[[length(count_draws_list) + 1]] <- tibble(
    Strategy = strat_name,
    Transition = "No ICD -> SCD",
    Value = sim_noicd_scd
  )
  count_draws_list[[length(count_draws_list) + 1]] <- tibble(
    Strategy = strat_name,
    Transition = "No ICD -> Death",
    Value = sim_noicd_death
  )
}

count_draws <- bind_rows(count_draws_list)

# 5. Compute posterior predictive p-values (ppp-value)
# ppp = P(Y_rep >= Y_obs)
ppp_stats <- count_draws %>%
  inner_join(raw_counts, by = c("Strategy", "Transition")) %>%
  group_by(Strategy, Transition) %>%
  summarise(
    Observed_Count = first(Observed_Count),
    Mean_Predicted = mean(Value),
    Median_Predicted = median(Value),
    Lower_95_Predicted = quantile(Value, 0.025),
    Upper_95_Predicted = quantile(Value, 0.975),
    ppp = mean(Value >= Observed_Count),
    .groups = "drop"
  )

print("----- Posterior Predictive Check Statistics -----")
print(as.data.frame(ppp_stats))

# 6. Plot 1: Posterior Transition Probabilities vs. Raw proportions
# Use ggdist for beautiful distribution plots
p1 <- ggplot(prob_draws, aes(x = Value, y = Transition, fill = Transition)) +
  stat_slab(alpha = 0.6, color = "black", size = 0.2) +
  stat_pointinterval(position = position_dodge(width = 0.4), .width = c(0.66, 0.95)) +
  geom_point(data = raw_probs, aes(x = Observed_Val, y = Transition),
             color = "#d95f02", size = 4, shape = 18, inherit.aes = FALSE) +
  xlim(0, 0.035) +
  facet_wrap(~ Strategy, ncol = 2, scales = "free_x") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Posterior Transition Probabilities vs. Raw Proportions",
    subtitle = "Posterior distributions (slab/interval) compared to raw proportions (orange diamonds).",
    x = "Transition Probability (per year)",
    y = "Transition Pathway",
    fill = "Transition"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 5)),
    plot.subtitle = element_text(face = "italic", size = 12, margin = margin(b = 15)),
    strip.background = element_rect(fill = "#2c3e50", color = NA),
    strip.text = element_text(face = "bold", color = "white", size = 12),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.5, "lines")
  )

p1

ggsave(filename = here::here("images/posterior_check_probabilities.png"), plot = p1,
       dpi = 600, width = 12, height = 9)

# 7. Plot 2: Posterior Predictive Checks (Counts)
# Create label data frame for facets
text_labels <- ppp_stats %>%
  mutate(
    label = paste0("Obs = ", Observed_Count, "\nppp = ", round(ppp, 2)),
    # Use Inf for positioning
    x = Inf,
    y = Inf
  )

p2 <- ggplot(count_draws, aes(x = Value, fill = Transition)) +
  geom_histogram(binwidth = 1, alpha = 0.7, color = "black", size = 0.2) +
  geom_vline(data = raw_counts, aes(xintercept = Observed_Count),
             color = "#e74c3c", linetype = "longdash", size = 1) +
  geom_text(data = text_labels, aes(x = x, y = y, label = label),
            hjust = 1.1, vjust = 1.5, size = 4.5, fontface = "bold",
            color = "#2c3e50", inherit.aes = FALSE) +
  facet_grid(Strategy ~ Transition, scales = "free") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Posterior Predictive Checks: Transition Counts",
    subtitle = "Simulated transition counts (histogram) compared to actual observed counts (red dashed line).",
    x = "Predicted Count of Events (per year)",
    y = "MCMC Iteration Frequency",
    fill = "Transition"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 5)),
    plot.subtitle = element_text(face = "italic", size = 12, margin = margin(b = 15)),
    strip.background = element_rect(fill = "#2c3e50", color = NA),
    strip.text = element_text(face = "bold", color = "white", size = 11),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines")
  )

p2

ggsave(filename = here::here("images/posterior_predictive_checks.png"), plot = p2,
       dpi = 600, width = 14, height = 10)
