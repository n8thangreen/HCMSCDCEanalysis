
#' convert BUGS transition probability matrix
#' to model with tunnel state
#'
#' @param lambda S by S matrix
#' @param n_tunnel Number of tunnel states (maximum device lifespan in cycles, default: 15)
#' @param alpha proportion move from hcm to icd
#' @param p_repl Numeric vector of length \code{n_tunnel}: probability of
#'   device replacement in each tunnel state (time since implant), e.g. from
#'   \code{\link{calc_repl_probs}}. The last element is forced to 1 so that
#'   all devices are replaced by the maximum lifespan. Default \code{NULL} gives
#'   deterministic replacement after \code{n_tunnel} cycles.
#'
lambda_tunnel <- function(lambda, n_tunnel = 15, alpha = 0, p_repl = NULL) {
  # icd and shock state are duplicated
  n_states <- n_tunnel*2 + 4
  trans_mat <- matrix(0, ncol = n_states, nrow = n_states)

  if (is.null(p_repl)) {
    p_repl <- c(rep(0, n_tunnel - 1), 1)
  }
  if (length(p_repl) != n_tunnel) {
    stop("`p_repl` must have length `n_tunnel` (", n_tunnel, ").")
  }
  if (anyNA(p_repl) || any(p_repl < 0 | p_repl > 1)) {
    stop("`p_repl` must contain probabilities in [0, 1].")
  }
  p_repl[n_tunnel] <- 1

  pdeath <- lambda[1,3]

  for (t in seq_len(n_tunnel)) {
    # next tunnel state; replacement resets the clock to state 1
    t_next <- if (t < n_tunnel) t + 1 else 1

    # icd -> icd: replaced (state 1) or ageing (state t + 1)
    trans_mat[t, 1] <- trans_mat[t, 1] + lambda[1,1] * p_repl[t]
    trans_mat[t, t_next] <- trans_mat[t, t_next] + lambda[1,1] * (1 - p_repl[t])

    # icd -> shock (same time since implant)
    trans_mat[t, n_tunnel + t] <- lambda[1,2]

    # shock -> icd: replaced or ageing
    trans_mat[n_tunnel + t, 1] <-
      trans_mat[n_tunnel + t, 1] + (1 - pdeath) * p_repl[t]
    trans_mat[n_tunnel + t, t_next] <-
      trans_mat[n_tunnel + t, t_next] + (1 - pdeath) * (1 - p_repl[t])
  }

  # all-cause mortality
  trans_mat[1:(2*n_tunnel), (2*n_tunnel) + 1] <- pdeath

  # sink states
  trans_mat[(2*n_tunnel) + 1, 2*n_tunnel + 1] <- 1
  trans_mat[(2*n_tunnel) + 3, 2*n_tunnel + 3] <- 1
  trans_mat[(2*n_tunnel) + 4, 2*n_tunnel + 4] <- 1

  # hcm
  trans_mat[(2*n_tunnel) + 2, 1] <- alpha*lambda[4,4]
  trans_mat[(2*n_tunnel) + 2, 2*n_tunnel + 2] <- (1 - alpha)*lambda[4,4]
  trans_mat[(2*n_tunnel) + 2, 2*n_tunnel + 3] <- lambda[4,5]
  trans_mat[(2*n_tunnel) + 2, 2*n_tunnel + 4] <- lambda[4,6]

  trans_mat
}
