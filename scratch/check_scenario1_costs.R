library(dplyr)
library(purrr)
library(BCEA)
library(readr)

devtools::load_all(here::here())

load("data/jags_output_all_strat_direct.RData")
lambda_samples <- mm1$BUGSoutput$sims.list$lambda

load("data/transdat.RData")
data("ipd_risk")

interv_names <- names(transdat)
n_intervs <- length(interv_names)
n_sim <- mm1$BUGSoutput$n.sims
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

inp_tab <- readr::read_csv(
  here::here("inst/extdata/scenario_input_values.csv"),
  col_names = TRUE, show_col_types = FALSE
)

# Scenario 1 (tvicd, u_icd = 0.9, q_hcm = 0.88)
inp <- inp_tab[1, ]
inp$t_repl_max <- 1

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

n_sample <- 3672
cols_to_compare <- c("risk_over_4", "risk_over_6", "icd")
costs <- colMeans(out$c_total[, cols_to_compare] / n_sample)
qalys <- colMeans(out$e_total[, cols_to_compare] / n_sample)

df_res <- data.frame(
  Strategy = cols_to_compare,
  Cost = costs,
  QALYs = qalys
)

cat("=== SCENARIO 1 (BASE-CASE) DETERMINISTIC RESULTS ===\n")
print(df_res)
