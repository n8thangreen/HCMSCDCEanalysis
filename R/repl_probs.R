#' Weibull scale parameter from a median
#'
#' @param median Median device lifespan (years)
#' @param shape Weibull shape parameter
#' @return Weibull scale parameter
#' @export
weibull_scale <- function(median, shape) {
  median / (log(2)^(1 / shape))
}

#' Cycle-specific device replacement probabilities from a Weibull curve
#'
#' Discrete-time hazard \eqn{1 - S(t)/S(t-1)} for each tunnel state
#' (time since implant), where \eqn{S} is the Weibull battery survival function.
#' The last element is set to 1 (forced replacement at maximum lifespan).
#'
#' @param max_cycles Number of tunnel states (\code{t_repl_max})
#' @param shape Weibull shape parameter
#' @param scale Weibull scale parameter
#' @return Numeric vector of length \code{max_cycles}
#' @export
calc_repl_probs <- function(max_cycles, shape, scale) {
  cycles <- seq_len(max_cycles)
  surv_t <- stats::pweibull(cycles, shape = shape, scale = scale, lower.tail = FALSE)
  surv_t_minus_1 <- stats::pweibull(cycles - 1, shape = shape, scale = scale, lower.tail = FALSE)

  p <- ifelse(surv_t_minus_1 == 0, 1, 1 - surv_t / surv_t_minus_1)
  p[max_cycles] <- 1
  p
}

#' Sample replacement probabilities for PSA
#'
#' Draws the Weibull median (Normal, truncated above 0) and shape (Gamma) for
#' each simulation and converts to a matrix of per-tunnel-state replacement
#' probabilities, analogous to the posterior samples of the other transition
#' probabilities. With zero (or \code{NULL}) standard errors every row is
#' identical (deterministic curve).
#'
#' @param n_sim Number of simulations
#' @param median Mean median device lifespan (years)
#' @param shape Mean Weibull shape
#' @param t_max Number of tunnel states (\code{t_repl_max})
#' @param median_se Standard error of the median
#' @param shape_se Standard error of the shape
#' @return Matrix \code{n_sim} x \code{t_max} of replacement probabilities
#' @export
sample_repl_probs <- function(n_sim, median, shape, t_max,
                              median_se = NULL, shape_se = NULL) {
  has_se <- function(se) !is.null(se) && !is.na(se) && se > 0

  med <- if (has_se(median_se)) {
    pmax(stats::rnorm(n_sim, median, median_se), 0.1)
  } else {
    rep(median, n_sim)
  }
  shp <- if (has_se(shape_se)) {
    stats::rgamma(n_sim, shape = (shape / shape_se)^2, rate = shape / shape_se^2)
  } else {
    rep(shape, n_sim)
  }

  t(vapply(seq_len(n_sim), function(i) {
    calc_repl_probs(t_max, shp[i], weibull_scale(med[i], shp[i]))
  }, numeric(t_max)))
}

