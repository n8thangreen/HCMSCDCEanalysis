test_that("weibull_scale gives S(median) = 0.5", {
  sc <- weibull_scale(median = 8.5, shape = 5)
  expect_equal(stats::pweibull(8.5, shape = 5, scale = sc, lower.tail = FALSE), 0.5)
})

test_that("calc_repl_probs returns hazards with forced replacement in last state", {
  sc <- weibull_scale(10, 5)
  p <- calc_repl_probs(15, shape = 5, scale = sc)
  expect_length(p, 15)
  expect_true(all(p >= 0 & p <= 1))
  expect_equal(p[15], 1)
  # hazard increases with device age (shape > 1)
  expect_true(all(diff(p[1:14]) > 0))
  # cumulative survival at the median is 0.5
  expect_equal(prod(1 - p[1:10]), 0.5)
})

test_that("sample_repl_probs has correct shape and deterministic limit", {
  m <- sample_repl_probs(20, median = 10, shape = 5, t_max = 15,
                         median_se = 0.5, shape_se = 0.5)
  expect_equal(dim(m), c(20, 15))
  expect_true(all(m[, 15] == 1))
  expect_gt(max(apply(m, 2, sd)[1:14]), 0)

  d <- sample_repl_probs(5, median = 10, shape = 5, t_max = 15)
  expect_true(all(apply(d, 2, sd) == 0))
  expect_equal(d[1, ], calc_repl_probs(15, 5, weibull_scale(10, 5)))
})

test_that("lambda_tunnel with p_repl builds hazard-based transitions", {
  lambda <- matrix(0, 6, 6)
  lambda[1, ] <- c(0.90, 0.05, 0.05, 0, 0, 0)
  lambda[4, ] <- c(0, 0, 0, 0.9, 0.05, 0.05)
  n_tunnel <- 4
  p_repl <- c(0.1, 0.2, 0.5, 1)

  m <- lambda_tunnel(lambda, n_tunnel = n_tunnel, alpha = 0, p_repl = p_repl)

  expect_equal(rowSums(m), rep(1, 2 * n_tunnel + 4))
  # ICD_t -> replace (state 1) or age (t + 1)
  expect_equal(m[2, 1], 0.90 * 0.2)
  expect_equal(m[2, 3], 0.90 * 0.8)
  # Shock_t -> replace or age
  expect_equal(m[n_tunnel + 3, 1], 0.95 * 0.5)
  expect_equal(m[n_tunnel + 3, 4], 0.95 * 0.5)
  # last state always replaced
  expect_equal(m[4, 1], 0.90)
  expect_equal(m[2 * n_tunnel, 1], 0.95)

  # default reproduces deterministic replacement
  m_det <- lambda_tunnel(lambda, n_tunnel = n_tunnel)
  expect_equal(m_det, lambda_tunnel(lambda, n_tunnel = n_tunnel, p_repl = c(0, 0, 0, 1)))
  expect_equal(m_det[1, 2], 0.90)
  expect_equal(m_det[1, 1], 0)

  expect_error(lambda_tunnel(lambda, n_tunnel = 4, p_repl = c(0.1, 0.2)), "length")
  expect_error(lambda_tunnel(lambda, n_tunnel = 4, p_repl = c(0.1, 0.2, 1.5, 1)), "probabilities")
})

test_that("create_probs_array accepts per-simulation p_repl", {
  n_sim <- 3
  lam <- array(0, dim = c(n_sim, 6, 6))
  lam[, 1, ] <- rep(c(0.9, 0.05, 0.05, 0, 0, 0), each = n_sim)
  lam[, 4, ] <- rep(c(0, 0, 0, 0.9, 0.05, 0.05), each = n_sim)
  lam_list <- list(a = lam, b = lam)
  p <- sample_repl_probs(n_sim, 3, 5, 5, median_se = 0.3, shape_se = 0.5)

  probs <- create_probs_array(lam_list, n_tunnel = 5, p_repl = p)
  expect_equal(dim(probs$a), c(14, 14, 1, n_sim))
  expect_equal(probs$a[1, 1, 1, 2], 0.9 * p[2, 1])
  expect_equal(probs$b[1, 2, 1, 3], 0.9 * (1 - p[3, 1]))

  probs_l <- create_probs_array(lam_list, n_tunnel = 5, p_repl = list(a = p, b = p * 0 + c(0, 0, 0, 0, 1)))
  expect_equal(probs_l$b[1, 1, 1, 1], 0)

  expect_error(create_probs_array(lam_list, n_tunnel = 5, p_repl = p[1:2, ]), "n_sim")
})

test_that("get_t_repl_max falls back to legacy t_repl with a warning", {
  expect_equal(get_t_repl_max(list(t_repl_max = 7)), 7)
  expect_warning(v <- get_t_repl_max(list(t_repl = 6)), "deprecated")
  expect_equal(v, 6)
})
