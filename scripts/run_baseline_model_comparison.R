#

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

# Find n_sim from jags_output
temp_env <- new.env()
load(here::here("data/jags_output.RData"), envir = temp_env)
n_sim <- temp_env$mm1$BUGSoutput$n.sims

J <- 30 # Time horizon
WTP_THRESHOLD <- 30000

# Read scenario parameter inputs
inp_tab <- readr::read_csv(
  here::here("inst/extdata/scenario_input_values.csv"),
  col_names = TRUE, show_col_types = FALSE
)

run_simulation <- function(device, model_type) {
  # Load BUGS output locally in function environment
  load(here::here("data/jags_output.RData"))

  lambda0_local <- mm1$BUGSoutput$sims.list$lambda.0
  lambda1_local <- mm1$BUGSoutput$sims.list$lambda.1
  lambda2_local <- mm1$BUGSoutput$sims.list$lambda.2
  lambda3_local <- mm1$BUGSoutput$sims.list$lambda.3

  # 1. Determine scenario label and filter
  lbl <- tolower(device)
  inp <- inp_tab %>%
    dplyr::filter(Scenario_Label == lbl) %>%
    dplyr::slice(1)

  # 2. Setup initial counts
  init_counts <- get_init_counts(ipd_risk)

  # 3. Build model components
  devtools::load_all(here::here())

  if (model_type == "static") {
    S <- 6
    inp$t_repl <- 1

    # Create starting population arrays
    pop <- create_pop_array(init_counts, J, n_sim, model_type = "static", interv_names = interv_names)

    # Create transition probability arrays
    lambda_list <- list(obs = lambda0_local, risk6 = lambda1_local, risk4 = lambda2_local, no_icd = lambda3_local)
    probs <- create_probs_array(lambda_list, model_type = "static")

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
  } else {
    n_tunnel <- inp$t_repl
    S <- 2 * n_tunnel + 4

    # Create starting population arrays
    pop <- create_pop_array(init_counts, J, n_sim, model_type = "tunnel", n_tunnel = n_tunnel, interv_names = interv_names)

    # Create transition probability arrays
    alpha <- 0
    lambda_list <- list(obs = lambda0_local, risk6 = lambda1_local, risk4 = lambda2_local, no_icd = lambda3_local)
    probs <- create_probs_array(lambda_list, model_type = "tunnel", n_tunnel = n_tunnel, alpha = alpha)

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
  }

  # Calculate mean QALYs and costs directly from simulation results (keep positive)
  qalys <- out$e_total / n_total
  costs <- out$c_total / n_total

  mean_qaly <- colMeans(qalys)
  mean_cost <- colMeans(costs)

  # Reference strategy is "No ICD" (named "no_icd")
  ref_qaly <- mean_qaly["no_icd"]
  ref_cost <- mean_cost["no_icd"]

  inc_qaly <- mean_qaly - ref_qaly
  inc_cost <- mean_cost - ref_cost
  inb <- (inc_qaly * WTP_THRESHOLD) - inc_cost

  # Extract statistics
  res_summary <- data.frame(
    Device = device,
    Model = model_type,
    Strategy = c("Risk > 4%", "Risk > 6%", "Observed", "No ICD"),
    Mean_QALY = c(mean_qaly["risk4"], mean_qaly["risk6"], mean_qaly["obs"], mean_qaly["no_icd"]),
    Mean_Cost = c(mean_cost["risk4"], mean_cost["risk6"], mean_cost["obs"], mean_cost["no_icd"]),
    Incremental_QALY = c(inc_qaly["risk4"], inc_qaly["risk6"], inc_qaly["obs"], inc_qaly["no_icd"]),
    Incremental_Cost = c(inc_cost["risk4"], inc_cost["risk6"], inc_cost["obs"], inc_cost["no_icd"]),
    INB = c(inb["risk4"], inb["risk6"], inb["obs"], inb["no_icd"])
  )
  return(res_summary)
}

# Run all 4 combinations
res_tv_std <- run_simulation("TVICD", "static")
res_tv_tun <- run_simulation("TVICD", "tunnel")
res_sicd_std <- run_simulation("SICD", "static")
res_sicd_tun <- run_simulation("SICD", "tunnel")

all_results <- bind_rows(res_tv_std, res_tv_tun, res_sicd_std, res_sicd_tun)

print(all_results, row.names = FALSE)

# Save to CSV
write_csv(all_results, here::here("data/all_ce_model_results.csv"))
