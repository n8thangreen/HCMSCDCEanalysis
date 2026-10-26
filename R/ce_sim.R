
#' Markov model cost-effectiveness simulation
#'
#' @param pop Number of individuals. List-matrix \code{n_init} states x time x n_sim
#'   Simulation parameters are taken from characteristics of pop array.
#' @param probs Transition probabilities. List-array \code{n_init} states x states x time x sim
#' @param inp Input parameters dataframe/list
#' @param c_unit Unit costs per state. List \code{t_max}\code{n_init}
#' @param e_unit Health unit values per state. List \code{n_init}
#' @param c_init One-off initial state costs. List \code{n_init}
#' @param e_init One-off initial state health values. List \code{n_init}
#' @param distn Logical. Sample from parameter distributions?
#' @param interv_names Names of interventions
#' @param delta Discount rate
#' @param ... Additional arguments passed to methods
#' @return List
#'
#' @importFrom purrr map
#'
#' @export
#'
ce_sim <- function(pop,
                   probs,
                   inp = NULL,
                   c_unit,
                   e_unit,
                   c_init = NA,
                   e_init = NA,
                   distn = TRUE,
                   interv_names = NULL,
                   delta = 0.035,
                   ...) {

  S <- dim(pop[[1]])[1]        # number of states
  tmax <- dim(pop[[1]])[2]     # time horizon
  n_sim <- dim(pop[[1]])[3]
  n_interv <- length(pop)
  n_tunnel <- if (S >= 6 && (S - 4) %% 2 == 0) (S - 4) / 2 else NULL
  if (is.null(interv_names)) {
    dots <- list(...)
    if ("labels" %in% names(dots)) {
      interv_names <- dots$labels
    } else {
      interv_names <- seq_len(n_interv)
    }
  }

  # initialise empty output matrices
  out_mat <-
    map(1:n_interv,
        ~matrix(NA,
                nrow = n_sim,
                ncol = tmax)) %>%
    stats::setNames(interv_names)

  cost <- out_mat
  eff <- out_mat
  dcost <- out_mat
  deff <- out_mat

  if (is.null(inp)) {
    distn <- FALSE
    n_params <- 0
  } else {
    n_params <- ifelse(distn, nrow(inp), ncol(inp))
  }
  input_tab <- matrix(NA, nrow = n_sim, ncol = n_params)

  for (i in seq_len(n_sim)) {

    # sample inputs
    input_vals <-
      if (distn) {sample_param(inp)} else {inp}

    if (n_params > 0) {
      input_tab[i, ] <- as.matrix(input_vals)
    }

    current_delta <- if (!is.null(input_vals) && !is.null(input_vals$delta)) {
      input_vals$delta
    } else {
      delta
    }

    for (k in interv_names) {

      # calculate state unit values
      rc_init <- safe_call(c_init, k, tmax, input_vals, n_tunnel = n_tunnel)
      re_init <- safe_call(e_init, k, tmax, input_vals, n_tunnel = n_tunnel)
      rc_unit <- safe_call(c_unit, k, tmax, input_vals, n_tunnel = n_tunnel)
      re_unit <- safe_call(e_unit, k, tmax, input_vals, n_tunnel = n_tunnel)

      for (j in seq_len(tmax)) {

        if (j > 1) {
          for (s in seq_len(S)) {
            pop[[k]][s, j, i] <- state_next(pop, probs, i, j, k, s)
          }
        }

        disc <- (1 + current_delta)^(j-1)

        cost[[k]][i, j] <- rc_unit[, j] %*% pop[[k]][, j, i]
        eff[[k]][i, j] <- re_unit[, j] %*% pop[[k]][, j, i]

        dcost[[k]][i, j] <- cost[[k]][i, j] / disc
        deff[[k]][i, j] <- eff[[k]][i, j] / disc
      }

      # add one-off starting state values
      if (is.function(c_init)) {
        c_init_tot <- rc_init %*% pop[[k]][, 1, i]

        cost[[k]][i, 1] <- cost[[k]][i, 1] + c_init_tot
        dcost[[k]][i, 1] <- dcost[[k]][i, 1] + c_init_tot
      }
      if (is.function(e_init)) {
        e_init_tot <- re_init %*% pop[[k]][, 1, i]

        eff[[k]][i, 1] <- eff[[k]][i, 1] + e_init_tot
        deff[[k]][i, 1] <- deff[[k]][i, 1] + e_init_tot
      }
    }
  }

  # total cost across all time points
  c_total <-
    purrr::map_dfc(dcost, rowSums) %>%
    as.matrix()

  e_total <-
    purrr::map_dfc(deff, rowSums) %>%
    as.matrix()

  list(pop = pop,
       cost = cost,
       c_total = c_total,
       dcost = dcost,
       eff = eff,
       deff = deff,
       e_total = e_total,
       input_tab = input_tab)
}


#
state_next <- function(pop, probs, i, j, k, s, t) {

  # time-homogeneous model has dim 1
  t <- min(dim(probs[[1]])[3], j - 1)

  t(pop[[k]][, j - 1, i]) %*% probs[[k]][, s, t, i]
}


# Helper to safely call a function with only the arguments it accepts
safe_call <- function(fn, k = NULL, ..., n_tunnel = NULL) {
  if (!is.function(fn)) return(NA)
  args <- list(...)
  fn_formals <- formals(fn)
  is_closure_zero_args <- !is.primitive(fn) && is.null(fn_formals)
  res <- if (is_closure_zero_args) {
    do.call(fn, list())
  } else if (is.null(fn_formals) || "..." %in% names(fn_formals)) {
    do.call(fn, args)
  } else {
    fn_arg_names <- names(fn_formals)
    call_args <- list()
    if ("J" %in% fn_arg_names) call_args$J <- args[[1]]
    if ("inp" %in% fn_arg_names) {
      call_args$inp <- if (length(args) >= 2) args[[2]] else args[[1]]
    }
    if ("name_interv" %in% fn_arg_names) call_args$name_interv <- k
    if ("n_tunnel" %in% fn_arg_names && !is.null(n_tunnel)) {
      call_args$n_tunnel <- n_tunnel
    }

    if (length(call_args) == 0 && length(fn_arg_names) > 0) {
      n_args <- length(fn_formals)
      call_args <- args[seq_len(min(length(args), n_args))]
    }
    do.call(fn, call_args)
  }

  if (is.list(res) && !is.null(k)) {
    if (!is.null(names(res)) && as.character(k) %in% names(res)) {
      res[[as.character(k)]]
    } else {
      k_idx <- suppressWarnings(as.numeric(k))
      if (!is.na(k_idx) && k_idx <= length(res)) {
        res[[k_idx]]
      } else {
        res
      }
    }
  } else {
    res
  }
}
