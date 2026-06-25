# Script to replicate paper results (Table 2)
# Using TV-ICD (device="ICD"), standard 6-state model, 12-year time horizon, and "observed" as reference.

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
J <- 12 # 12-year time horizon as in the paper

# 3. Read baseline scenario parameters
inp_tab <- readr::read_csv(
  here::here("inst/extdata/scenario_input_values.csv"),
  col_names = TRUE, show_col_types = FALSE
)

# Filter for TV-ICD (scenario 1 in scenario_input_values.csv)
inp <- inp_tab %>%
  dplyr::filter(Scenario_Label == "tvicd") %>%
  dplyr::slice(1)

# Enforce no tunnel states for standard model
inp$t_repl <- 1

# Setup initial starting populations from individual patient data
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

# Standard 6-state model starting populations
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

# Copy background mortality (ICD -> Death) to Shock state (Shock -> Death)
lambda.0[, STATE_SHOCK, STATE_DEATH] <- lambda.0[, STATE_ICD, STATE_DEATH]
lambda.1[, STATE_SHOCK, STATE_DEATH] <- lambda.1[, STATE_ICD, STATE_DEATH]
lambda.2[, STATE_SHOCK, STATE_DEATH] <- lambda.2[, STATE_ICD, STATE_DEATH]
lambda.3[, STATE_SHOCK, STATE_DEATH] <- lambda.3[, STATE_ICD, STATE_DEATH]

# Set Shock -> ICD transition as the remainder
lambda.0[, STATE_SHOCK, STATE_ICD]   <- 1 - lambda.0[, STATE_SHOCK, STATE_DEATH]
lambda.1[, STATE_SHOCK, STATE_ICD]   <- 1 - lambda.1[, STATE_SHOCK, STATE_DEATH]
lambda.2[, STATE_SHOCK, STATE_ICD]   <- 1 - lambda.2[, STATE_SHOCK, STATE_DEATH]
lambda.3[, STATE_SHOCK, STATE_ICD]   <- 1 - lambda.3[, STATE_SHOCK, STATE_DEATH]

# Probability transition lists
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

# Calculate BCEA results using 'obs' as the reference
he <- bcea(
  eff = -out$e_total / n_total,
  cost = -out$c_total / n_total,
  ref = 3, # Reference is 'obs' (3rd element in interv_names)
  interventions = c("Cox 5-year SCD risk >4%",
                    "Cox 5-year SCD risk >6%",
                    "Observed Strategy",
                    "No ICD Strategy")
)

# Print summary
cat("\n=== REPLICATED PAPER TABLE 2 SUMMARY ===\n")
print(summary(he))

# Print mean cost and QALY for each strategy
cat("\n=== ABSOLUTE VALUES ===\n")
abs_costs <- colMeans(out$c_total / n_total)
abs_qalys <- colMeans(out$e_total / n_total)
names(abs_costs) <- c("Risk > 4%", "Risk > 6%", "Observed", "No ICD")
names(abs_qalys) <- c("Risk > 4%", "Risk > 6%", "Observed", "No ICD")

df_abs <- data.frame(
  Strategy = names(abs_costs),
  Cost = abs_costs,
  QALYs = abs_qalys
)
print(df_abs)
