# Baseline Model Comparison Script (TV-ICD vs. S-ICD Tunnel-State Model)
#
# Runs the tunnel-state Markov model for TV-ICD (10-year replacement cycle)
# and S-ICD (8-year replacement cycle) over a 30-year horizon and outputs
# summary statistics to data/all_ce_model_results.csv.

library(dplyr)
library(purrr)
library(BCEA)
library(readr)
library(here)

devtools::load_all(here::here())

# Load data
load(here::here("data/transdat.RData"))
data("ipd_risk")

n_total <- nrow(ipd_risk)
interv_names <- c("risk4", "risk6", "obs", "no_icd")

# Load BUGS output
load(here::here("data/jags_output.RData"))
lambda0_local <- mm1$BUGSoutput$sims.list$lambda.0
lambda1_local <- mm1$BUGSoutput$sims.list$lambda.1
lambda2_local <- mm1$BUGSoutput$sims.list$lambda.2
lambda3_local <- mm1$BUGSoutput$sims.list$lambda.3
n_sim <- mm1$BUGSoutput$n.sims

J <- 30 # Time horizon (years)
WTP_THRESHOLD <- 25000

# Read scenario parameter inputs
inp_tab <- readr::read_csv(
  here::here("inst/extdata/scenario_input_values.csv"),
  col_names = TRUE, show_col_types = FALSE
)

run_simulation <- function(device) {
  lbl <- tolower(device)
  inp <- inp_tab %>%
    dplyr::filter(Scenario_Label == lbl) %>%
    dplyr::slice(1)

  init_counts <- get_init_counts(ipd_risk)
  p_repl <- sample_repl_probs(n_sim, inp$repl_median, inp$repl_shape)  # deterministic Weibull (no SEs)

  pop <- create_pop_array(init_counts, J, n_sim, interv_names = interv_names)

  lambda_list <- list(obs = lambda0_local, risk6 = lambda1_local, risk4 = lambda2_local, no_icd = lambda3_local)
  probs <- create_probs_array(lambda_list, alpha = 0, p_repl = p_repl)

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

  qalys <- out$e_total / n_total
  costs <- out$c_total / n_total

  mean_qaly <- colMeans(qalys)
  mean_cost <- colMeans(costs)

  # Reference strategy within device is "Observed" ("obs") and also report vs "No ICD" ("no_icd")
  ref_qaly <- mean_qaly["obs"]
  ref_cost <- mean_cost["obs"]

  inc_qaly <- mean_qaly - ref_qaly
  inc_cost <- mean_cost - ref_cost
  icer <- ifelse(inc_qaly == 0, NA_real_, inc_cost / inc_qaly)
  inb <- (inc_qaly * WTP_THRESHOLD) - inc_cost

  res_summary <- data.frame(
    Device = device,
    Model = "tunnel (15yr)",
    Strategy = c("Risk > 4%", "Risk > 6%", "Observed", "No ICD"),
    Mean_QALY = c(mean_qaly["risk4"], mean_qaly["risk6"], mean_qaly["obs"], mean_qaly["no_icd"]),
    Mean_Cost = c(mean_cost["risk4"], mean_cost["risk6"], mean_cost["obs"], mean_cost["no_icd"]),
    Incremental_QALY = c(inc_qaly["risk4"], inc_qaly["risk6"], inc_qaly["obs"], inc_qaly["no_icd"]),
    Incremental_Cost = c(inc_cost["risk4"], inc_cost["risk6"], inc_cost["obs"], inc_cost["no_icd"]),
    ICER = c(icer["risk4"], icer["risk6"], icer["obs"], icer["no_icd"]),
    INB = c(inb["risk4"], inb["risk6"], inb["obs"], inb["no_icd"])
  )
  return(res_summary)
}

res_tv_tun <- run_simulation("TVICD")
res_sicd_tun <- run_simulation("SICD")

all_results <- bind_rows(res_tv_tun, res_sicd_tun)

print(all_results, row.names = FALSE)

write_csv(all_results, here::here("data/all_ce_model_results.csv"))
