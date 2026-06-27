# =============================================================================
# Script to run Deterministic Sensitivity Analysis (DSA)
# Generates data/res_dsa.RData for one-way tornado plots
# =============================================================================

library(dplyr)
library(purrr)
library(readr)
devtools::load_all(here::here())

# 1. Load Data & MCMC Output
# -----------------------------------------------------------------------------
data("transdat")
data("ipd_risk")

load(here::here("data/jags_output.RData"))

# Explicitly extract transition probabilities (avoiding attach.bugs())
trans_prob_obs    <- mm1$BUGSoutput$sims.list$lambda.0
trans_prob_risk6  <- mm1$BUGSoutput$sims.list$lambda.1
trans_prob_risk4  <- mm1$BUGSoutput$sims.list$lambda.2
trans_prob_no_icd <- mm1$BUGSoutput$sims.list$lambda.3

n_sim   <- mm1$BUGSoutput$n.sims
n_total <- nrow(ipd_risk)
J       <- 12 # 12-year time horizon for base-case
S       <- 6  # Standard memoryless model

interv_names <- c("risk4", "risk6", "obs", "no_icd")

# 2. Setup Initial Populations
# -----------------------------------------------------------------------------
init <- inits_by_interv(ipd_risk, group_names = c("icd", "risk_over_6", "risk_over_4"))

init_counts <- list(
  obs    = list(ICD = dplyr::filter(init, risk == "icd", status == "ICD")$n,
                no_ICD = dplyr::filter(init, risk == "icd", status == "low_risk")$n),
  risk6  = list(ICD = dplyr::filter(init, risk == "risk_over_6", status == "ICD")$n,
                no_ICD = dplyr::filter(init, risk == "risk_over_6", status == "low_risk")$n),
  risk4  = list(ICD = dplyr::filter(init, risk == "risk_over_4", status == "ICD")$n,
                no_ICD = dplyr::filter(init, risk == "risk_over_4", status == "low_risk")$n),
  no_icd = list(ICD = 0,
                no_ICD = n_total)
)

n_init <- list(
  obs    = c(init_counts$obs$ICD,    0, 0, init_counts$obs$no_ICD,    0, 0),
  risk6  = c(init_counts$risk6$ICD,  0, 0, init_counts$risk6$no_ICD,  0, 0),
  risk4  = c(init_counts$risk4$ICD,  0, 0, init_counts$risk4$no_ICD,  0, 0),
  no_icd = c(init_counts$no_icd$ICD, 0, 0, init_counts$no_icd$no_ICD, 0, 0)
)

# 3. Format Transition Probabilities
# -----------------------------------------------------------------------------

# Map probabilities into the standard 6-state arrays
probs_empty <- array(NA, dim = c(S, S, 1, n_sim))
probs <- list(obs = probs_empty, risk6 = probs_empty, risk4 = probs_empty, no_icd = probs_empty)

probs$obs[,, 1, ]    <- aperm(trans_prob_obs, c(2, 3, 1))
probs$risk6[,, 1, ]  <- aperm(trans_prob_risk6, c(2, 3, 1))
probs$risk4[,, 1, ]  <- aperm(trans_prob_risk4, c(2, 3, 1))
probs$no_icd[,, 1, ] <- aperm(trans_prob_no_icd, c(2, 3, 1))

# 4. Run the DSA Loop
# -----------------------------------------------------------------------------
inp_tab <- readr::read_csv(
  here::here("inst/extdata/scenario_input_values.csv"),
  show_col_types = FALSE
)

res <- list()

for (i in 1:nrow(inp_tab)) {
  inp <- inp_tab[i, ]

  # Enforce standard 6-state model behavior (no tunnel states)
  inp$t_repl <- 1

  # Initialize fresh populations for this scenario
  pop <- list()
  for (k in interv_names) {
    pop[[k]] <- array(NA, dim = c(S, J, n_sim))
    pop[[k]][, 1, ] <- t(n_init[[k]])
  }

  # Run deterministic simulation (distn = FALSE)
  res[[i]] <- ce_sim(
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
}

# 5. Save Output
save(res, file = here::here("data/res_dsa.RData"))
