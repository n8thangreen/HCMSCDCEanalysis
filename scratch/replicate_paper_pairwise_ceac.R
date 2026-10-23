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
S <- 3

p_death_bg <- lambda_samples[,, 1, 3]
lambda_samples[,, 2, 3] <- p_death_bg
lambda_samples[,, 2, 1] <- 1 - p_death_bg

probs <- list()
for (j in seq_along(interv_names)) {
  nm <- interv_names[j]
  probs_array <- array(NA, dim = c(2 * S, 2 * S, 1, n_sim))
  interv_slice <- lambda_samples[, j, , ]
  probs_array[,, 1, ] <- aperm(interv_slice, c(2, 3, 1))
  probs[[nm]] <- probs_array
}

# Setup starting populations
init_counts <- inits_by_interv(ipd_risk, group_names = interv_names) |>
  mutate(state = "hcm")

init_full <- expand.grid(
  risk = interv_names,
  status = c("ICD", "low_risk"),
  state = c("hcm", "scd/shock", "death")) |>
  arrange(risk, status) |>
  left_join(init_counts, by = c("risk", "status", "state")) |>
  tidyr::replace_na(list(n = 0))

pop0 <- split(init_full$n, init_full$risk)
pop <- list()
for (i in interv_names) {
  pop[[i]] <- array(NA, dim = c(2*S, J, n_sim))
  pop[[i]][, 1, ] <- pop0[[i]]
}

# 3. Read input_distns.csv and filter for vticd
inp_dist <- readr::read_csv(
  here::here("inst/extdata/input_distns.csv"),
  col_names = TRUE, show_col_types = FALSE
)

inp <- inp_dist %>%
  dplyr::filter(label == "vticd")

# Enforce paper's 2 appointments, 0.043 complication rate, no tunnel states, zero replacement cost, and fix health parameters
inp <- inp %>%
  mutate(params = case_when(
    variable == "c_appt" ~ "290",       # Point at 290
    variable == "p_compl" ~ "0.043",    # Point at 0.043
    variable == "t_repl" ~ "1",         # Set t_repl to 1
    variable == "c_icd_repl" ~ "0",     # Set replacement cost to 0
    variable == "q_hcm" ~ "0.88",       # Fix health parameters
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

# Run simulation
out <- ce_sim(
  pop = pop,
  probs = probs,
  inp = inp,
  c_unit = rc_unit,
  e_unit = re_unit,
  c_init = rc_init,
  e_init = re_init,
  distn = TRUE, # PSA
  interv_names = interv_names
)

# Extract only the 3 strategies compared in the paper: risk_over_4, risk_over_6, and icd
cols_to_compare <- c("risk_over_4", "risk_over_6", "icd")
eff_subset <- out$e_total[, cols_to_compare]
cost_subset <- out$c_total[, cols_to_compare]

# Reference is "icd" (Observed Strategy), which is 3rd column in the subset
# We run pairwise comparisons against Observed (icd)
he_6 <- bcea(
  eff = eff_subset[, c("risk_over_6", "icd")] / n_total,
  cost = cost_subset[, c("risk_over_6", "icd")] / n_total,
  ref = 2,
  interventions = c("Cox 5-year SCD risk >6%", "Observed Strategy")
)

he_4 <- bcea(
  eff = eff_subset[, c("risk_over_4", "icd")] / n_total,
  cost = cost_subset[, c("risk_over_4", "icd")] / n_total,
  ref = 2,
  interventions = c("Cox 5-year SCD risk >4%", "Observed Strategy")
)

# Extract CEAC probabilities
k_vals <- he_6$k
ceac_6 <- he_6$ceac[, 1] # Probability that risk_over_6 is cost-effective vs Observed
ceac_4 <- he_4$ceac[, 1] # Probability that risk_over_4 is cost-effective vs Observed

df_ceac <- data.frame(
  k = rep(k_vals, 2),
  Probability = c(ceac_6, ceac_4),
  Strategy = rep(c("Cox 5-year SCD risk >6%", "Cox 5-year SCD risk >4%"), each = length(k_vals))
)

# Plot CEAC
gg_ceac <- ggplot(df_ceac, aes(x = k, y = Probability, linetype = Strategy)) +
  geom_line(linewidth = 1, color = "black") +
  scale_linetype_manual(values = c("Cox 5-year SCD risk >6%" = "solid", "Cox 5-year SCD risk >4%" = "dashed")) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey55") +
  geom_vline(xintercept = 25000, linetype = "dashed", color = "grey55") +
  labs(
    title = "Cost-Effectiveness Acceptability Curve",
    x = "Willingness-to-pay threshold k (£)",
    y = "Probability of cost-effectiveness"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

# Save plot
ggsave(
  filename = "C:/Users/Nathan/.gemini/antigravity-cli/brain/25b9f518-fdb8-43cd-8051-d98663134580/ceac_replicated_pairwise.png",
  plot = gg_ceac,
  width = 6,
  height = 6,
  dpi = 320
)


# Save CE Plane using BCEA on the subset of 3 strategies (Observed is reference)
# Since they are plotted relative to Observed, the standard plot shows the incremental cost & QALYs
he_both <- bcea(
  eff = -eff_subset / n_total,
  cost = -cost_subset / n_total,
  ref = 3, # Observed is 3rd column
  interventions = c("Cox 5-year SCD risk >4%",
                    "Cox 5-year SCD risk >6%",
                    "Observed Strategy")
)

ceac.plot(he_both, graph = "ggplot")

png(filename = "C:/Users/Nathan/.gemini/antigravity-cli/brain/25b9f518-fdb8-43cd-8051-d98663134580/ce_plane_replicated_pairwise.png", width = 640, height = 640, res = 120)
ceplane.plot(he_both, graph = "base", title = "Cost-Effectiveness Plane")
dev.off()

cat("Pairwise plots generated successfully.\n")
