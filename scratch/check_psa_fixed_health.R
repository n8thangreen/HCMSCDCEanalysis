library(dplyr)
library(purrr)
library(BCEA)
library(readr)

devtools::load_all(here::here())

# 1. Load transition counts and patient risk profiles
data("transdat")
data("ipd_risk")

# 2. Load the BUGS MCMC output
load(here::here("data/jags_output.RData"))
R2WinBUGS::attach.bugs(mm1$BUGSoutput) 

n_total <- nrow(ipd_risk)
interv_names <- c("risk4", "risk6", "obs", "no_icd")
n_sim <- mm1$BUGSoutput$n.sims
J <- 12 # 12-year time horizon

# Setup initial starting populations
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

# Set Shock transitions
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

pop <- list()
for (i in interv_names) {
  pop[[i]] <- array(NA, dim = c(S, J, n_sim))
  pop[[i]][, 1, ] <- t(n_init[[i]])
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
    variable == "t_repl_max" ~ "1",         # Set t_repl_max to 1
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

cols_to_compare <- c("risk4", "risk6", "obs")
costs <- colMeans(out$c_total[, cols_to_compare] / n_total)
qalys <- colMeans(out$e_total[, cols_to_compare] / n_total)

df_res <- data.frame(
  Strategy = cols_to_compare,
  Cost = costs,
  QALYs = qalys
)

cat("=== PSA RESULTS WITH FIXED HEALTH PARAMETERS ===\n")
print(df_res)
