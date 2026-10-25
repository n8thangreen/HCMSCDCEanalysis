# Replicate paper using jags_output_all_strat_direct.RData transition probabilities
library(dplyr)
library(purrr)
library(BCEA)
library(readr)
library(ggplot2)

devtools::load_all(here::here())

# 1. Load transition counts and patient risk profiles
load(here::here("data/transdat.RData"))
data("ipd_risk")

# 2. Load the BUGS MCMC output
load(here::here("data/jags_output_all_strat_direct.RData"))
lambda_samples <- mm1$BUGSoutput$sims.list$lambda
n_sim <- dim(lambda_samples)[1]
interv_names <- names(transdat)
n_intervs <- length(interv_names)

n_total <- nrow(ipd_risk)
J <- 12 # 12-year time horizon
S <- 3  # state dimension in BUGS (hcm, scd/shock, death)

# Process probabilities to 2*S = 6 states total
probs <- list()
p_death_bg <- lambda_samples[,, 1, 3]
lambda_samples[,, 2, 3] <- p_death_bg
lambda_samples[,, 2, 1] <- 1 - p_death_bg

for (j in seq_along(interv_names)) {
  nm <- interv_names[j]
  probs_array <- array(NA, dim = c(2 * S, 2 * S, 1, n_sim))
  interv_slice <- lambda_samples[, j, , ]
  probs_array[,, 1, ] <- aperm(interv_slice, c(2, 3, 1))
  probs[[nm]] <- probs_array
}

# Setup initial populations
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

# Load input_distns.csv and filter for vticd
inp_dist <- readr::read_csv(
  here::here("inst/extdata/input_distns.csv"),
  col_names = TRUE, show_col_types = FALSE
)

inp <- inp_dist %>%
  dplyr::filter(label == "vticd")

# Enforce paper's parameters
inp <- inp %>%
  mutate(params = case_when(
    variable == "c_appt" ~ "290",       # Point at 290
    variable == "p_compl" ~ "0.043",    # Point at 0.043
    variable == "t_repl_max" ~ "1",         # Set t_repl_max to 1 (enforcing no tunnel states)
    variable == "c_icd_repl" ~ "0",     # Set replacement cost to 0
    TRUE ~ params
  ))

rpoint <- function(n, x) x

set.seed(42)

# Run the standard 6-state simulation with distn = TRUE
out <- ce_sim(
  pop = pop,
  probs = probs,
  inp = inp,
  c_unit = rc_unit,
  e_unit = re_unit,
  c_init = rc_init,
  e_init = re_init,
  distn = TRUE,
  interv_names = interv_names
)

# Extract only the 3 strategies compared in the paper's Figure 2: risk_over_4, risk_over_6, and icd
# Order of columns in out$e_total corresponds to interv_names:
# "risk_over_6", "risk_over_4", "any_rf", "over_1_rf", "over_2_rf", "over_3_rf", "icd", "no_icd"
cols_to_compare <- c("risk_over_4", "risk_over_6", "icd")
eff_subset <- out$e_total[, cols_to_compare]
cost_subset <- out$c_total[, cols_to_compare]

# Reference is "icd" (Observed Strategy), which is 3rd column in the subset
he <- bcea(
  eff = -eff_subset / n_total,
  cost = -cost_subset / n_total,
  ref = 3,
  interventions = c("Cox 5-year SCD risk >4%",
                    "Cox 5-year SCD risk >6%",
                    "Observed Strategy")
)

cat("\n========================================================\n")
cat("=== PSA RESULTS WITH DIRECT JAGS OUTPUT ===\n")
cat("========================================================\n")
print(summary(he))

abs_costs <- colMeans(cost_subset / n_total)
abs_qalys <- colMeans(eff_subset / n_total)
df_abs <- data.frame(
  Strategy = cols_to_compare,
  Cost = abs_costs,
  QALYs = abs_qalys
)
print(df_abs)

# 2. Generate and save Cost-effectiveness Plane
cat("\nGenerating Cost-effectiveness Plane plot (base R)...\n")
png(filename = here::here("images/ce_plane_original_icd_replicated_direct_jags.png"), width = 640, height = 640, res = 120)
ceplane.plot(he, graph = "base", title = "Cost-Effectiveness Plane")
dev.off()

# 3. Generate and save Cost-effectiveness Acceptability Curve (CEAC)
cat("Generating Cost-effectiveness Acceptability Curve (CEAC) plot (base R)...\n")
png(filename = here::here("images/ceac_original_icd_replicated_direct_jags.png"), width = 640, height = 640, res = 120)
ceac.plot(he, graph = "base", title = "Cost-Effectiveness Acceptability Curve")
dev.off()

cat("\nPlots saved successfully to images/.\n")
