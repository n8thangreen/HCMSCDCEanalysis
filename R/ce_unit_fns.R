
#' Sample unit costs (per cycle) with tunnel states
#'
#' @param J Maximum time horizon (year), numeric
#' @param inp Input parameters
#' @param name_interv Name of intervention
#'
#' @returns matrix
#' @export
rc_unit <- function(J, inp, name_interv) {

  n_tunnel <- get_t_repl_max(inp)

  # --- 1. ICD ARM (Tunnel States) ---
  # Standard Appt cost applies to all years
  v_icd <- rep(inp$c_appt, n_tunnel)

  # --- 2. EVENTS & OTHER STATES ---
  v_shock <- rep(inp$c_shock, n_tunnel) # Duplicated to match shock tunnel states
  v_death <- 0
  v_scd   <- 0

  # --- 3. NO-ICD ARM ---
  # Usually zero cost for healthy No-ICD patients (single state)
  v_noicd <- 0

  # --- COMBINE ---
  # Order must match the transition matrix:
  # [ICD_Tunnel (T), Shock_Tunnel (T), Death_ICD (1), NoICD_Healthy (1), SCD (1), Death_NoICD (1)]

  cost_vector <- c(
    v_icd,      # ICD Tunnel (size T)
    v_shock,    # Shock Tunnel (size T)
    v_death,    # Death ICD (size 1)
    v_noicd,    # No-ICD Healthy (size 1)
    v_scd,      # SCD State (size 1)
    v_death     # Death No-ICD (size 1)
  )

  # (Rows = States, Cols = Time)
  cost_mat <- matrix(cost_vector, nrow = length(cost_vector), ncol = J, byrow = FALSE)

  # ADD LATE LEAD COMPLICATION COSTS (yearly, by lead age = model cycle)
  # Expected cost of late infection and dislodgement/failure while holding an ICD
  # (ICD and shock tunnel states). Not applied to the no-ICD strategy.
  if (name_interv != "no_icd") {
    idx_dev <- seq_len(2 * n_tunnel)
    lead_cost <- get_param(inp, "p_inf_yr") * get_param(inp, "c_inf") +
      p_lead_dis(inp, J) * get_param(inp, "c_dis")
    cost_mat[idx_dev, ] <- sweep(cost_mat[idx_dev, , drop = FALSE], 2, lead_cost, "+")
  }

  # ADD REPLACEMENT COST
  # Applied to the FIRST tunnel state only for cycles j > 1 (when n_tunnel > 1).
  # At j = 1, patients are newly implanted and their initial device cost (c_icd)
  # is accounted for in rc_init(). When the device expires at the end of state T,
  # the transition matrix sends surviving patients back to State 1 in cycle j > 1.
  # Replacement carries its own complication risk (p_repl_compl).
  if (n_tunnel > 1 && J > 1) {
    cost_mat[1, 2:J] <- cost_mat[1, 2:J] + inp$c_icd_repl +
      get_param(inp, "p_repl_compl") * c_compl_weighted(inp)
  }

  return(cost_mat)
}

#' Get an optional input parameter
#'
#' @param inp Input parameters
#' @param name Parameter name
#' @param default Value returned if the parameter is absent
#' @returns numeric scalar
#' @keywords internal
get_param <- function(inp, name, default = 0) {
  val <- inp[[name]]
  if (is.null(val) || length(val) == 0 || is.na(val)) default else as.numeric(val)
}

#' Get the maximum device lifespan (number of tunnel states)
#'
#' Returns \code{inp$t_repl_max}. For backwards compatibility falls back to the
#' legacy name \code{inp$t_repl} (with a warning).
#'
#' @param inp Input parameters
#' @returns numeric scalar
#' @export
get_t_repl_max <- function(inp) {
  val <- inp[["t_repl_max"]]
  if (is.null(val) && !is.null(inp[["t_repl"]])) {
    warning("`t_repl` is deprecated; use `t_repl_max`.", call. = FALSE)
    val <- inp[["t_repl"]]
  }
  as.numeric(val)
}
#' Weighted cost of a lead complication (infection or dislodgement)
#'
#' @param inp Input parameters
#' @returns numeric scalar; 0 if both initial probabilities are 0
#' @export
c_compl_weighted <- function(inp) {
  p_inf <- get_param(inp, "p_inf_init")
  p_dis <- get_param(inp, "p_dis_init")
  denom <- p_inf + p_dis
  if (denom > 0) {
    (p_inf * get_param(inp, "c_inf") + p_dis * get_param(inp, "c_dis")) / denom
  } else {
    0
  }
}

#' Yearly probability of late lead dislodgement/failure by lead age
#'
#' Weibull cumulative hazard H(a) = lambda * a^gamma with lambda set from the
#' year-1 probability \code{p_dis_yr}. Lead age equals the model cycle because
#' all ICDs are implanted at cycle 1 and leads are not changed at generator
#' replacement. \code{gamma_dis = 1} gives a constant yearly probability.
#' The hazard stops rising after \code{t_lead_cap} years of lead age (it stays at
#' the year-\code{t_lead_cap} value), because the Weibull is not supported by
#' data beyond that horizon. Default \code{Inf} means no cap.
#'
#' @param inp Input parameters (\code{p_dis_yr}, \code{gamma_dis}, \code{t_lead_cap})
#' @param J Number of cycles
#' @returns numeric vector of length J
#' @export
p_lead_dis <- function(inp, J) {
  p1 <- get_param(inp, "p_dis_yr")
  g <- get_param(inp, "gamma_dis", default = 1)
  cap <- get_param(inp, "t_lead_cap", default = Inf)
  if (p1 <= 0) return(numeric(J))
  lambda <- -log(1 - p1)
  a <- seq_len(J)
  a_eff <- pmin(a, max(1, floor(cap)))
  1 - exp(-lambda * (a_eff^g - (a_eff - 1)^g))
}


#' Sample starting costs (at t=0)
#' @param inp Input parameters
#' @param name_interv Name of intervention
#'
#' different for each intervention scenario
#' only != 0 for hcm icd state
#' @return numeric vector
#' @export
rc_init <- function(inp, name_interv) {

  n_tunnel <- get_t_repl_max(inp)

  # Initialize empty vector of correct length
  # (2 * n_tunnel + 4)
  v_init <- numeric(2 * n_tunnel + 4)

  # --- CALCULATE COMPONENT COSTS ---

  # 1. Implant Procedure (Weighted by complication probability)
  # If prob is 0, cost is 0 to avoid division by zero error
  c_compl_w <- c_compl_weighted(inp)

  # Total immediate cost of getting an ICD
  c_procedure <- inp$c_icd + (inp$p_compl * c_compl_w)

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
  idx_noicd_start <- 2 * n_tunnel + 2  # T (ICD) + T (Shock) + 1 (Death) + 1 (Start NoICD)

  if (name_interv == "no_icd") {
    # Strategy: No ICD.
    # Cost: Only Testing (c_rscore).
    # Patient starts in: No-ICD Healthy state
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
#' @return matrix
#' @export
re_unit <- function(J, inp, name_interv) {

  n_tunnel <- get_t_repl_max(inp)

  # Utilities
  u_hcm   <- inp$q_hcm
  u_icd   <- u_hcm * inp$u_icd   # Utility of living with ICD
  u_shock <- u_icd * inp$u_shock # Utility decrement for Shock

  # --- BUILD VECTOR ---

  # ICD Arm
  v_icd   <- rep(u_icd, n_tunnel)
  v_shock <- rep(u_shock, n_tunnel) # Duplicated to match shock tunnel states

  # No-ICD Arm
  v_noicd <- u_hcm # Single state

  # Full Vector
  # [ICD_Tunnel (T), Shock_Tunnel (T), Death_ICD (1), NoICD_Healthy (1), SCD (1), Death_NoICD (1)]

  util_vector <- c(
    v_icd,     # ICD Healthy
    v_shock,   # Shock Tunnel
    0,         # Death ICD
    v_noicd,   # No-ICD Healthy
    0,         # SCD (Assumed 0 QALY)
    0          # Death No-ICD
  )

  util_mat <- matrix(util_vector,
                     nrow = length(util_vector),
                     ncol = J,
                     byrow = FALSE)

  # ADD LATE LEAD COMPLICATION DISUTILITY (yearly, by lead age = model cycle)
  # Expected disutility of late infection and dislodgement/failure while holding
  # an ICD (ICD and shock tunnel states). Not applied to the no-ICD strategy.
  if (name_interv != "no_icd") {
    idx_dev <- seq_len(2 * n_tunnel)
    lead_disutil <- (get_param(inp, "p_inf_yr") + p_lead_dis(inp, J)) *
      get_param(inp, "u_compl")
    util_mat[idx_dev, ] <- sweep(util_mat[idx_dev, , drop = FALSE], 2, lead_disutil, "+")
  }

  # ADD REPLACEMENT DISUTILITY
  # Patients re-entering tunnel state 1 after device expiry (cycles j > 1) receive
  # a new ICD, so incur the implant decrement with the replacement complication
  # probability (see re_init()). Mirrors the replacement cost logic in rc_unit().
  if (name_interv != "no_icd" && n_tunnel > 1 && J > 1) {
    util_mat[1, 2:J] <- util_mat[1, 2:J] + implant_disutility(inp, replacement = TRUE)
  }

  return(util_mat)
}

#' Disutility incurred on receiving a new ICD (initial implant or replacement)
#'
#' @param inp Input parameters
#' @param replacement Logical. If TRUE use the replacement decrement
#'   \code{u_repl} (defaults to \code{u_implant} if absent) and the replacement
#'   complication probability \code{p_repl_compl} instead of \code{u_implant}
#'   and \code{p_compl}.
#' @returns numeric scalar
#' @export
implant_disutility <- function(inp, replacement = FALSE) {
  if (replacement) {
    get_param(inp, "u_repl", default = inp$u_implant) +
      get_param(inp, "p_repl_compl") * inp$u_compl
  } else {
    inp$u_implant + inp$p_compl * inp$u_compl
  }
}

#' Sample starting disutilities (at t=0)
#' @param inp Input parameters
#' @param name_interv Name of intervention
#' @returns numeric vector
#' @export
re_init <- function(inp, name_interv) {

  n_tunnel <- get_t_repl_max(inp)
  v_init <- numeric(2 * n_tunnel + 4)

  # Calculate Disutility
  # u_implant is usually negative
  e_decrement <- implant_disutility(inp)

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


#' Quantile function for a point-mass (fixed value) "distribution"
#'
#' @param p Probabilities (ignored, only the length is used)
#' @param x The fixed value
#' @returns numeric vector the same length as \code{p}
#' @export
qpoint <- function(p, x) rep(x, length(p))

#' sample from generic distribution
#'
#' Draws one uniform per parameter and transforms it with the distribution's
#' quantile function (\code{q<distn>}). The marginal distributions are the same
#' as with \code{r<distn>}, but a fixed number of random numbers is consumed
#' whatever the distribution family, so runs for different device labels with
#' the same seed share common random numbers (paired comparisons are not
#' swamped by Monte Carlo noise).
#'
#' @param inp Input parameters
#' @export
sample_param <- function(inp) {
  u <- stats::runif(length(inp$variable))
  sample_fns <- paste0("q", inp$distn, "(", sprintf("%.17g", u), ", ", inp$params, ")")
  rparam <- map(sample_fns, ~eval(parse(text = .)))
  names(rparam) <- inp$variable
  as.data.frame(rparam)
}
