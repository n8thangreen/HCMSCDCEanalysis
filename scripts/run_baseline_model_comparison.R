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

  # 3. Build model components
  devtools::load_all(here::here())

  if (model_type == "static") {
    S <- 6
    inp$t_repl <- 1

    n_init <- list(
      obs    = c(init_counts$obs$ICD,    0, 0, init_counts$obs$no_ICD,    0, 0),
      risk6  = c(init_counts$risk6$ICD,  0, 0, init_counts$risk6$no_ICD,  0, 0),
      risk4  = c(init_counts$risk4$ICD,  0, 0, init_counts$risk4$no_ICD,  0, 0),
      no_icd = c(init_counts$no_icd$ICD, 0, 0, init_counts$no_icd$no_ICD, 0, 0)
    )

    probs_empty <- array(NA, dim = c(S, S, 1, n_sim))
    probs <- list(obs = probs_empty, risk6 = probs_empty, risk4 = probs_empty, no_icd = probs_empty)

    probs$obs[,, 1, ]    <- aperm(lambda0_local, c(2, 3, 1))
    probs$risk6[,, 1, ]  <- aperm(lambda1_local, c(2, 3, 1))
    probs$risk4[,, 1, ]  <- aperm(lambda2_local, c(2, 3, 1))
    probs$no_icd[,, 1, ] <- aperm(lambda3_local, c(2, 3, 1))

    pop <- list()
    for (i in interv_names) {
      pop[[i]] <- array(NA, dim = c(S, J, n_sim))
      pop[[i]][, 1, ] <- t(n_init[[i]])
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
  } else {
    n_tunnel <- inp$t_repl
    S <- 2 * n_tunnel + 4

    probs_empty <- array(NA, dim = c(S, S, 1, n_sim))
    probs <- list(obs = probs_empty, risk6 = probs_empty, risk4 = probs_empty, no_icd = probs_empty)

    alpha <- 0
    for (i in 1:n_sim) {
      probs$obs[,, 1, i]    <- lambda_tunnel(lambda0_local[i,,], n_tunnel = n_tunnel, alpha)
      probs$risk6[,, 1, i]  <- lambda_tunnel(lambda1_local[i,,], n_tunnel = n_tunnel, alpha)
      probs$risk4[,, 1, i]  <- lambda_tunnel(lambda2_local[i,,], n_tunnel = n_tunnel, alpha)
      probs$no_icd[,, 1, i] <- lambda_tunnel(lambda3_local[i,,], n_tunnel = n_tunnel, alpha)
    }

    pop <- list()
    for (i in interv_names) {
      pop[[i]] <- array(NA, dim = c(S, J, n_sim))
      pop[[i]][, 1, ] <- 0
      pop[[i]][1, 1, ] <- init_counts[[i]]$ICD
      pop[[i]][2 * n_tunnel + 2, 1, ] <- init_counts[[i]]$no_ICD
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
