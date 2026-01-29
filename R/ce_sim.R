
#' Markov model cost-effectiveness simulation
#'
#' @param pop Number of individuals. List-matrix \code{n_init} states x time x n_sim
#'   Simulation parameters are taken from characteristics of pop array.
#' @param probs Transition probabilities. List-array \code{n_init} states x states x time x sim
#' @param c_unit Unit costs per state. List \code{t_max}\code{n_init}
#' @param e_unit Health unit values per state. List \code{n_init}
#' @param c_init One-off initial state costs. List \code{n_init}
#' @param interv_names Names of interventions
#' @return List
#' @seealso \code{\link{init_pop}}
#'
#' @importFrom purrr map
#'
#' @export
#'
ce_sim <- function(pop,
                   probs,
                   inp,
                   c_unit,
                   e_unit,
                   c_init = NA,
                   e_init = NA,
                   distn = TRUE,
                   interv_names = NULL) {

  S <- dim(pop[[1]])[1]        # number of states
  tmax <- dim(pop[[1]])[2]     # time horizon
  n_sim <- dim(pop[[1]])[3]
  n_interv <- length(pop)
  if (is.null(interv_names)) interv_names <- seq_len(n_interv)

  # initialise empty output matrices
  out_mat <-
    map(1:n_interv,
        ~matrix(NA,
                nrow = n_sim,
                ncol = tmax)) %>%
    setNames(interv_names)

  cost <- out_mat
  eff <- out_mat
  dcost <- out_mat
  deff <- out_mat

  n_params <- ifelse(distn, nrow(inp), ncol(inp))
  input_tab <- matrix(NA, nrow = n_sim, ncol = n_params)

  for (i in seq_len(n_sim)) {

    # sample inputs
    input_vals <-
      if (distn) {sample_param(inp)} else {inp}

    input_tab[i, ] <- as.matrix(input_vals)

    delta <- input_vals$delta

    for (k in interv_names) {

      # calculate state unit values
      rc_init <- c_init(input_vals, k)
      re_init <- e_init(input_vals, k)
      rc_unit <- c_unit(tmax, input_vals, k)
      re_unit <- e_unit(tmax, input_vals, k)

      for (j in seq_len(tmax)) {

        if (j > 1) {
          for (s in seq_len(S)) {
            pop[[k]][s, j, i] <- state_next(pop, probs, i, j, k, s)
          }
        }

        disc <- (1 + delta)^(j-1)

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

