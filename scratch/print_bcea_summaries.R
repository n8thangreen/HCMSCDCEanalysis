library(dplyr)
library(purrr)
library(BCEA)
library(readr)

devtools::load_all(here::here())

load(here::here("data/transdat.RData"))
data("ipd_risk")
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

inp_dist <- readr::read_csv(
  here::here("inst/extdata/input_distns.csv"),
  col_names = TRUE, show_col_types = FALSE
)

inp <- inp_dist %>%
  dplyr::filter(label == "vticd")

inp <- inp %>%
  mutate(params = case_when(
    variable == "c_appt" ~ "290",
    variable == "p_compl" ~ "0.043",
    variable == "t_repl" ~ "1",
    variable == "c_icd_repl" ~ "0",
    variable == "q_hcm" ~ "0.88",
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

cols_to_compare <- c("risk_over_4", "risk_over_6", "icd")
eff_subset <- out$e_total[, cols_to_compare]
cost_subset <- out$c_total[, cols_to_compare]

he_pos <- bcea(
  eff = eff_subset / n_total,
  cost = cost_subset / n_total,
  ref = 3,
  interventions = c("Cox 5-year SCD risk >4%",
                    "Cox 5-year SCD risk >6%",
                    "Observed Strategy")
)

he_neg <- bcea(
  eff = -eff_subset / n_total,
  cost = -cost_subset / n_total,
  ref = 3,
  interventions = c("Cox 5-year SCD risk >4%",
                    "Cox 5-year SCD risk >6%",
                    "Observed Strategy")
)

cat("=== BCEA SUMMARY (POSITIVE INPUTS) ===\n")
print(summary(he_pos))

cat("\n=== BCEA SUMMARY (NEGATIVE INPUTS) ===\n")
print(summary(he_neg))
