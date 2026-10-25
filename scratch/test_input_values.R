# Script to test parameters from input_values.csv
library(dplyr)
library(purrr)
library(BCEA)
library(readr)

devtools::load_all(here::here())

data("transdat")
data("ipd_risk")

load(here::here("data/jags_output.RData"))
R2WinBUGS::attach.bugs(mm1$BUGSoutput) 

n_total <- nrow(ipd_risk)
interv_names <- c("risk4", "risk6", "obs", "no_icd")
n_sim <- mm1$BUGSoutput$n.sims
J <- 12 # 12-year time horizon

init <- inits_by_interv(ipd_risk, group_names = c("icd", "risk_over_6", "risk_over_4"))
init_counts <- list(
  obs = list(
    ICD    = dplyr::filter(init, risk == "icd", status == "ICD")$n,
    no_ICD = dplyr::filter(init, risk == "icd", status == "low_risk")$n
  ),
  risk6 = list(
    ICD    = dplyr::filter(init, risk == "risk_over_6", status == "ICD")$n,
    no_ICD = dplyr::filter(init, risk == "risk_over_6", status == "low_risk")$n
  ),
  risk4 = list(
    ICD    = dplyr::filter(init, risk == "risk_over_4", status == "ICD")$n,
    no_ICD = dplyr::filter(init, risk == "risk_over_4", status == "low_risk")$n
  ),
  no_icd = list(
    ICD    = 0,
    no_ICD = n_total
  )
)

S <- 6
n_init <- list(
  obs    = c(init_counts$obs$ICD,    0, 0, init_counts$obs$no_ICD,    0, 0),
  risk6  = c(init_counts$risk6$ICD,  0, 0, init_counts$risk6$no_ICD,  0, 0),
  risk4  = c(init_counts$risk4$ICD,  0, 0, init_counts$risk4$no_ICD,  0, 0),
  no_icd = c(init_counts$no_icd$ICD, 0, 0, init_counts$no_icd$no_ICD, 0, 0)
)

STATE_ICD   <- 1
STATE_SHOCK <- 2
STATE_DEATH <- 3

lambda.0[, STATE_SHOCK, STATE_DEATH] <- lambda.0[, STATE_ICD, STATE_DEATH]
lambda.1[, STATE_SHOCK, STATE_DEATH] <- lambda.1[, STATE_ICD, STATE_DEATH]
lambda.2[, STATE_SHOCK, STATE_DEATH] <- lambda.2[, STATE_ICD, STATE_DEATH]
lambda.3[, STATE_SHOCK, STATE_DEATH] <- lambda.3[, STATE_ICD, STATE_DEATH]

lambda.0[, STATE_SHOCK, STATE_ICD]   <- 1 - lambda.0[, STATE_SHOCK, STATE_DEATH]
lambda.1[, STATE_SHOCK, STATE_ICD]   <- 1 - lambda.1[, STATE_SHOCK, STATE_DEATH]
lambda.2[, STATE_SHOCK, STATE_ICD]   <- 1 - lambda.2[, STATE_SHOCK, STATE_DEATH]
lambda.3[, STATE_SHOCK, STATE_ICD]   <- 1 - lambda.3[, STATE_SHOCK, STATE_DEATH]

probs_empty <- array(NA, dim = c(S, S, 1, n_sim))
probs <- list(obs = probs_empty, risk6 = probs_empty, risk4 = probs_empty, no_icd = probs_empty)

probs$obs[,, 1, ]    <- aperm(lambda.0, c(2, 3, 1))
probs$risk6[,, 1, ]  <- aperm(lambda.1, c(2, 3, 1))
probs$risk4[,, 1, ]  <- aperm(lambda.2, c(2, 3, 1))
probs$no_icd[,, 1, ] <- aperm(lambda.3, c(2, 3, 1))

# Read input_values.csv
inp_raw <- readr::read_csv(
  here::here("inst/extdata/input_values.csv"),
  col_names = TRUE, show_col_types = FALSE
)

# Convert to wide format (like scenario_input_values.csv)
inp <- inp_raw %>%
  dplyr::filter(label == "vticd") %>%
  dplyr::select(variable, value) %>%
  tidyr::pivot_wider(names_from = variable, values_from = value)

# Add other columns if they are missing
if (!"q_death" %in% names(inp)) inp$q_death <- 0
if (!"u_compl" %in% names(inp)) inp$u_compl <- -0.096
if (!"c_compl" %in% names(inp)) inp$c_compl <- 28857
if (!"delta" %in% names(inp)) inp$delta <- 0.035

# Enforce no tunnel states for standard model
inp$t_repl_max <- 1

pop <- list()
for (k in interv_names) {
  pop[[k]] <- array(NA, dim = c(S, J, n_sim))
  pop[[k]][, 1, ] <- t(n_init[[k]])
}

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

costs <- colMeans(out$c_total / n_total)
qalys <- colMeans(out$e_total / n_total)

df_res <- tibble(
  Strategy = names(costs),
  Cost = costs,
  QALYs = qalys
)

cat("=== RESULTS USING input_values.csv (vticd) ===\n")
print(df_res)
