library(dplyr)
library(purrr)
library(BCEA)
library(readr)

devtools::load_all(here::here())

# 1. Load transition counts and patient risk profiles
load(here::here("data/transdat.RData"))
data("ipd_risk")

# 2. Load the BUGS MCMC output
load(here::here("data/jags_output_all_strat_direct.RData"))
lambda_samples <- mm1$BUGSoutput$sims.list$lambda

# Get the mean transitions
lambda_mean <- apply(lambda_samples, c(2, 3, 4), mean)

n_total <- nrow(ipd_risk)
interv_names <- names(transdat)
n_sim <- 1
J <- 12 # 12-year time horizon
S <- 3  # state dimension in BUGS (hcm, scd/shock, death)

# Process probabilities to 2*S = 6 states total
probs <- list()
p_death_bg <- lambda_mean[, 1, 3] # ICD -> Death
lambda_mean[, 2, 3] <- p_death_bg
lambda_mean[, 2, 1] <- 1 - p_death_bg

for (j in seq_along(interv_names)) {
  nm <- interv_names[j]
  probs_array <- array(NA, dim = c(2 * S, 2 * S, 1, n_sim))
  probs_array[,, 1, 1] <- lambda_mean[j, , ]
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
  pop[[i]][, 1, 1] <- pop0[[i]]
}

# Load baseline scenario parameters
inp_tab <- readr::read_csv(
  here::here("inst/extdata/scenario_input_values.csv"),
  col_names = TRUE, show_col_types = FALSE
)

inp <- inp_tab %>%
  dplyr::filter(Scenario_Label == "tvicd") %>%
  dplyr::slice(1)

inp$t_repl <- 1
inp$c_appt <- 290
inp$p_compl <- 0.043

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

cols_to_compare <- c("risk_over_4", "risk_over_6", "icd")
costs <- colMeans(out$c_total[, cols_to_compare, drop=FALSE] / n_total)
qalys <- colMeans(out$e_total[, cols_to_compare, drop=FALSE] / n_total)

df_res <- tibble(
  Strategy = cols_to_compare,
  Cost = costs,
  QALYs = qalys
)

cat("=== DETERMINISTIC RESULTS WITH DIRECT JAGS OUTPUT ===\n")
print(df_res)
