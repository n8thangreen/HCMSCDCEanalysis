
#' Sample unit costs (per cycle) with tunnel states
#'
#' @param J Maximum time horizon (year), numeric
#' @param name_interv Name of intervention
#'
#' @returns list
#' @export
rc_unit <- function(J, inp, name_interv) {

  n_tunnel <- inp$t_repl

  # --- 1. ICD ARM (Tunnel States) ---
  # Standard Appt cost applies to all years
  v_icd <- rep(inp$c_appt, n_tunnel)

  # ADD REPLACEMENT COST
  # We apply it to the FIRST tunnel state.
  # Because when the device expires (State T), the transition matrix
  # sends the patient back to State 1. So State 1 represents "New Device Year".
  v_icd[1] <- v_icd[1] + inp$c_icd_repl

  # --- 2. EVENTS & OTHER STATES ---
  v_shock <- inp$c_shock # One-off cost for shock event
  v_death <- 0
  v_scd   <- 0

  # --- 3. NO-ICD ARM ---
  # Usually zero cost for healthy No-ICD patients
  v_noicd <- rep(0, n_tunnel)

  # --- COMBINE ---
  # Order must match the transition matrix:
  # [ICD_Tunnel (T), Shock (1), Death (1), NoICD_Tunnel (T), SCD (1), Death (1)]

  cost_vector <- c(
    v_icd,      # ICD Tunnel
    v_shock,    # Shock State
    v_death,    # Death ICD
    v_noicd,    # No-ICD Tunnel
    v_scd,      # SCD State
    v_death     # Death No-ICD
  )

  # (Rows = States, Cols = Time)
  # using replicate to fill across time J
  return(matrix(cost_vector, nrow = length(cost_vector), ncol = J, byrow = FALSE))
}


#' Sample starting costs (at t=0)
#' @param inp Input parameters
#' @param name_interv Name of intervention
#'
#' different for each intervention scenario
#' only != 0 for hcm icd state
#' @return list
#' @export
rc_init <- function(inp, name_interv) {

  n_tunnel <- inp$t_repl

  # Initialize empty vector of correct length
  # (2 * n_tunnel + 4)
  v_init <- numeric(2 * n_tunnel + 4)

  # --- CALCULATE COMPONENT COSTS ---

  # 1. Implant Procedure (Weighted by complication probability)
  # If prob is 0, cost is 0 to avoid division by zero error
  c_compl_weighted <- 0

  if ((inp$p_inf_init + inp$p_dis_init) > 0) {
    c_compl_weighted <- (inp$p_inf_init * inp$c_inf + inp$p_dis_init * inp$c_dis) /
      (inp$p_inf_init + inp$p_dis_init)
  }

  # Total immediate cost of getting an ICD
  c_procedure <- inp$c_icd + (inp$p_compl * c_compl_weighted)

  # 2. Risk Stratification Testing (c_rscore)
  # This applies if the strategy required testing patients
  is_testing_strategy <-
    name_interv %in% c("risk_over_4", "risk_over_6",
                       "any_rf", "over_1_rf", "over_2_rf", "over_3_rf",
                       "no_icd")

  c_test <- if (is_testing_strategy) inp$c_rscore else 0

  # --- ASSIGN TO CORRECT STATE ---

  # Define Start Indices
  idx_icd_start   <- 1
  idx_noicd_start <- n_tunnel + 3  # T (ICD) + 1 (Shock) + 1 (Death) + 1 (Start NoICD)

  if (name_interv == "no_icd") {
    # Strategy: No ICD.
    # Cost: Only Testing (c_rscore).
    # Patient starts in: No-ICD Tunnel State 1
    v_init[idx_noicd_start] <- c_test

  } else {
    # Strategy: ICD (Observed or Risk-based)
    # Cost: Implant + Testing (if applicable)
    # Patient starts in: ICD Tunnel State 1
    v_init[idx_icd_start] <- c_procedure + c_test
  }

  return(v_init)
}


#' Sample state utilities (QALYs per cycle) with tunnel states
#'
#' @param J Maximum time horizon (year)
#' @param inp Input parameters
#' @param name_interv Name of intervention
#' @return list
#' @export
re_unit <- function(J, inp, name_interv) {

  n_tunnel <- inp$t_repl

  # Utilities
  u_hcm   <- inp$q_hcm
  u_icd   <- u_hcm * inp$u_icd   # Utility of living with ICD
  u_shock <- u_icd * inp$u_shock # Utility decrement for Shock

  # --- BUILD VECTOR ---

  # ICD Arm
  v_icd   <- rep(u_icd, n_tunnel)

  # No-ICD Arm
  v_noicd <- rep(u_hcm, n_tunnel)

  # Full Vector
  # [ICD_Tunnel (T), Shock (1), Death (1), NoICD_Tunnel (T), SCD (1), Death (1)]

  util_vector <- c(
    v_icd,     # ICD Healthy
    u_shock,   # Shock State (Single state)
    0,         # Death
    v_noicd,   # No-ICD Healthy
    0,         # SCD (Assumed 0 QALY)
    0          # Death
  )

  return(matrix(util_vector,
                nrow = length(util_vector),
                ncol = J,
                byrow = FALSE))
}

#' Sample starting disutilities (at t=0)
#' @param inp Input parameters
#' @param name_interv Name of intervention
#' @returns list
#' @export
re_init <- function(inp, name_interv) {

  n_tunnel <- inp$t_repl
  v_init <- numeric(2 * n_tunnel + 4)

  # Calculate Disutility
  # u_implant is usually negative
  e_decrement <- inp$u_implant + (inp$p_compl * inp$u_compl)

  if (name_interv != "no_icd") {
    # Apply to ICD Start State (Index 1)
    v_init[1] <- e_decrement
  }

  return(v_init)
}

#' weighted infection, dislodgement
#' @param inp Input parameters
#' @export
c_complication <- function(inp) {

  ifelse(
    inp$p_inf_init + inp$p_dis_init == 0, 0,
    inp$p_inf_init*inp$c_inf + inp$p_dis_init*inp$c_dis)/(inp$p_inf_init + inp$p_dis_init)
}


#' sample from generic distribution
#' @param inp Input parameters
#' @export
sample_param <- function(inp) {
  sample_fns <- paste0("r", inp$distn, "(1, ", inp$params, ")")
  rparam <- map(sample_fns, ~eval(parse(text = .)))
  names(rparam) <- inp$variable
  as.data.frame(rparam)
}

