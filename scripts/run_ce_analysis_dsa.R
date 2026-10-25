# =============================================================================
# Script to run Deterministic Sensitivity Analysis (DSA) with Tunnel States
# Generates data/res_dsa.RData for one-way tornado plots (TV-ICD and S-ICD)
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

trans_prob_obs    <- mm1$BUGSoutput$sims.list$lambda.0
trans_prob_risk6  <- mm1$BUGSoutput$sims.list$lambda.1
trans_prob_risk4  <- mm1$BUGSoutput$sims.list$lambda.2
trans_prob_no_icd <- mm1$BUGSoutput$sims.list$lambda.3

n_sim   <- mm1$BUGSoutput$n.sims
n_total <- nrow(ipd_risk)
J       <- 30 # 30-year time horizon with replacement cycles

interv_names <- c("risk4", "risk6", "obs", "no_icd")
init_counts <- get_init_counts(ipd_risk)
lambda_list <- list(obs = trans_prob_obs, risk6 = trans_prob_risk6, risk4 = trans_prob_risk4, no_icd = trans_prob_no_icd)

# 2. Pre-build pop and probs by n_tunnel so we don't rebuild for each scenario row
# -----------------------------------------------------------------------------
inp_tab <- readr::read_csv(
  here::here("inst/extdata/scenario_input_values.csv"),
  show_col_types = FALSE
)

# cache key: tunnel length + Weibull replacement parameters (deterministic curve in DSA)
repl_key <- function(x) paste(x$t_repl_max, x$repl_median, x$repl_shape, sep = "_")
unique_tunnels <- unique(repl_key(inp_tab))
pop_cache <- list()
probs_cache <- list()

for (key in unique_tunnels) {
  row <- inp_tab[match(key, repl_key(inp_tab)), ]
  nt <- row$t_repl_max
  p_repl <- sample_repl_probs(n_sim, row$repl_median, row$repl_shape, nt)
  pop_cache[[key]] <- create_pop_array(init_counts, J, n_sim, n_tunnel = nt, interv_names = interv_names)
  probs_cache[[key]] <- create_probs_array(lambda_list, n_tunnel = nt, alpha = 0, p_repl = p_repl)
}

# 3. Run the DSA Loop across all 36 scenarios (1-18 TV-ICD, 19-36 S-ICD)
# -----------------------------------------------------------------------------
res <- list()

for (i in 1:nrow(inp_tab)) {
  inp <- inp_tab[i, ]
  key <- repl_key(inp)
  cat("Running DSA Scenario", i, "(", inp$Scenario_Label, ", t_repl_max =", inp$t_repl_max, ")...\n")

  res[[i]] <- ce_sim(
    pop = pop_cache[[key]],
    probs = probs_cache[[key]],
    inp = inp,
    c_unit = rc_unit,
    e_unit = re_unit,
    c_init = rc_init,
    e_init = re_init,
    distn = FALSE,
    interv_names = interv_names
  )
}

# 4. Save Output
save(res, file = here::here("data/res_dsa.RData"))
