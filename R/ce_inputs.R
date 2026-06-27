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
#' @param model_type Type of model: "static" or "tunnel".
#' @param n_tunnel Number of tunnel states (ignored if model_type is "static").
#' @param interv_names Names of the interventions (default: c("risk4", "risk6", "obs", "no_icd")).
#' @return A list of arrays, one per intervention, representing starting population states.
#' @export
create_pop_array <- function(init_counts, J, n_sim, model_type, n_tunnel = 1, interv_names = c("risk4", "risk6", "obs", "no_icd")) {
  S <- if (model_type == "static") 6 else (2 * n_tunnel + 4)
  pop <- list()
  for (i in interv_names) {
    pop[[i]] <- array(0, dim = c(S, J, n_sim))
    if (model_type == "static") {
      # For the static 6-state model, states are:
      # 1: ICD, 2: Shock, 3: Death ICD, 4: NoICD Healthy, 5: SCD, 6: Death NoICD
      pop[[i]][, 1, ] <- c(init_counts[[i]]$ICD, 0, 0, init_counts[[i]]$no_ICD, 0, 0)
    } else {
      # For the tunnel model, states are:
      # 1..T: ICD tunnel states, T+1..2T: Shock tunnel states, 2T+1: Death ICD, 2T+2: NoICD Healthy, 2T+3: SCD, 2T+4: Death NoICD
      pop[[i]][1, 1, ] <- init_counts[[i]]$ICD
      pop[[i]][2 * n_tunnel + 2, 1, ] <- init_counts[[i]]$no_ICD
    }
  }
  return(pop)
}

#' Create transition probability array from MCMC samples
#'
#' @param lambda_list List of transition probability arrays (one per intervention, size n_sim x 6 x 6).
#' @param model_type Type of model: "static" or "tunnel".
#' @param n_tunnel Number of tunnel states (ignored if model_type is "static").
#' @param alpha Proportion moving from hcm to icd (default: 0).
#' @return A list of transition probability arrays, one per intervention.
#' @export
create_probs_array <- function(lambda_list, model_type, n_tunnel = 1, alpha = 0) {
  n_sim <- dim(lambda_list[[1]])[1]
  S <- if (model_type == "static") 6 else (2 * n_tunnel + 4)
  probs_empty <- array(NA, dim = c(S, S, 1, n_sim))
  probs <- list()
  
  for (nm in names(lambda_list)) {
    probs[[nm]] <- probs_empty
    if (model_type == "static") {
      # Permute n_sim x 6 x 6 to 6 x 6 x 1 x n_sim
      probs[[nm]][,, 1, ] <- aperm(lambda_list[[nm]], c(2, 3, 1))
    } else {
      # For tunnel, loop over simulations and call lambda_tunnel
      for (i in 1:n_sim) {
        probs[[nm]][,, 1, i] <- lambda_tunnel(lambda_list[[nm]][i,,], n_tunnel = n_tunnel, alpha = alpha)
      }
    }
  }
  return(probs)
}
