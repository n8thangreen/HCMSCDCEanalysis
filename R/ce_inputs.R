#' Get initial state counts from individual patient data
#'
#' @param ipd_risk Dataframe containing patient risk data.
#' @return A list with counts of ICD and no_ICD patients for each strategy.
#' @importFrom dplyr filter
#' @export
get_init_counts <- function(ipd_risk) {
  init <- inits_by_interv(ipd_risk, group_names = c("icd", "risk_over_6", "risk_over_4"))
  n_total <- nrow(ipd_risk)
  list(
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
}

#' Create cohort starting population arrays
#'
#' @param init_counts List of initial state counts (ICD and no_ICD) per strategy.
#' @param J Time horizon (cycles/years).
#' @param n_sim Number of simulations.
#' @param n_tunnel Number of tunnel states (default: 1, which gives the 6-state non-tunnel model).
#' @param interv_names Names of the interventions (default: c("risk4", "risk6", "obs", "no_icd")).
#' @param ... Additional legacy arguments (e.g. `model_type`), ignored as non-tunnel is the `n_tunnel = 1` special case.
#' @return A list of arrays, one per intervention, representing starting population states.
#' @export
create_pop_array <- function(init_counts, J, n_sim, n_tunnel = 1, interv_names = c("risk4", "risk6", "obs", "no_icd"), ...) {
  S <- 2 * n_tunnel + 4
  pop <- list()

  for (i in interv_names) {
    pop[[i]] <- array(0, dim = c(S, J, n_sim))
    # States are:
    # 1..T: ICD tunnel states, T+1..2T: Shock tunnel states, 2T+1: Death ICD,
    # 2T+2: NoICD Healthy, 2T+3: SCD, 2T+4: Death NoICD
    pop[[i]][1, 1, ] <- init_counts[[i]]$ICD
    pop[[i]][2 * n_tunnel + 2, 1, ] <- init_counts[[i]]$no_ICD
  }
  return(pop)
}

#' Create transition probability array from MCMC samples
#'
#' @param lambda_list List of transition probability arrays (one per intervention, size n_sim x 6 x 6).
#' @param n_tunnel Number of tunnel states (default: 1, which gives the 6-state non-tunnel model).
#' @param alpha Proportion moving from hcm to icd (default: 0).
#' @param p_repl Device replacement probabilities by tunnel state: either an
#'   `n_sim` x `n_tunnel` matrix (e.g. from [sample_repl_probs()]) shared by all
#'   interventions, or a named list of such matrices (one per intervention).
#'   Default `NULL` gives deterministic replacement after `n_tunnel` cycles.
#' @param ... Additional legacy arguments (e.g. `model_type`), ignored as non-tunnel is the `n_tunnel = 1` special case.
#' @return A list of transition probability arrays, one per intervention.
#' @export
create_probs_array <- function(lambda_list, n_tunnel = 1, alpha = 0, p_repl = NULL, ...) {
  n_sim <- dim(lambda_list[[1]])[1]
  S <- 2 * n_tunnel + 4
  probs_empty <- array(NA, dim = c(S, S, 1, n_sim))
  probs <- list()

  for (nm in names(lambda_list)) {
    probs[[nm]] <- probs_empty
    p_repl_nm <- if (is.list(p_repl)) p_repl[[nm]] else p_repl
    if (!is.null(p_repl_nm) && nrow(p_repl_nm) < n_sim) {
      stop("`p_repl` must have at least `n_sim` (", n_sim, ") rows.")
    }
    for (i in 1:n_sim) {
      probs[[nm]][,, 1, i] <- lambda_tunnel(
        lambda_list[[nm]][i,,], n_tunnel = n_tunnel, alpha = alpha,
        p_repl = if (is.null(p_repl_nm)) NULL else p_repl_nm[i, ])
    }
  }
  return(probs)
}
