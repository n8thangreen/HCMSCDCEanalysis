inp <- data.frame(
  t_repl = 10, q_hcm = 0.88, u_icd = 0.9, u_shock = 0.875,
  u_implant = -0.048, u_compl = -0.096, u_repl = -0.024, p_compl = 0.043,
  p_repl_compl = 0.0091,
  c_appt = 145, c_shock = 1766, c_icd = 15317, c_icd_repl = 15317,
  c_inf = 37116, c_dis = 6146, p_inf_init = 0.02277, p_dis_init = 0.00828,
  p_inf_yr = 0, p_dis_yr = 0, gamma_dis = 1, c_rscore = 0
)

with_lead <- function(x, p_inf_yr = 0.002, p_dis_yr = 0.0029, gamma_dis = 1.6) {
  x$p_inf_yr <- p_inf_yr
  x$p_dis_yr <- p_dis_yr
  x$gamma_dis <- gamma_dis
  x
}

test_that("implant_disutility combines implant and complication decrements", {
  expect_equal(implant_disutility(inp), -0.048 + 0.043 * -0.096)
  expect_equal(implant_disutility(inp, replacement = TRUE), -0.024 + 0.0091 * -0.096)
  # falls back to u_implant if u_repl is absent
  x <- inp; x$u_repl <- NULL
  expect_equal(implant_disutility(x, replacement = TRUE), -0.048 + 0.0091 * -0.096)
})

test_that("re_unit applies disutility on replacement into tunnel state 1", {
  J <- 30
  mat <- re_unit(J, inp, "icd")
  base <- inp$q_hcm * inp$u_icd

  expect_equal(mat[1, 1], base)                       # initial implant handled by re_init
  expect_equal(mat[1, 2:J],
               rep(base + implant_disutility(inp, replacement = TRUE), J - 1))
  expect_equal(mat[2, ], rep(base, J))                # other tunnel states unchanged
})

test_that("re_unit applies no replacement disutility with a single tunnel state", {
  inp1 <- inp
  inp1$t_repl <- 1
  mat <- re_unit(10, inp1, "icd")
  expect_equal(mat[1, ], rep(inp1$q_hcm * inp1$u_icd, 10))
})

test_that("re_unit applies no replacement disutility for no_icd", {
  mat <- re_unit(10, inp, "no_icd")
  expect_equal(mat[1, ], rep(inp$q_hcm * inp$u_icd, 10))
})

test_that("p_lead_dis is constant when gamma_dis = 1 and increasing when > 1", {
  const <- p_lead_dis(with_lead(inp, gamma_dis = 1), 30)
  expect_equal(const, rep(0.0029, 30))

  inc <- p_lead_dis(with_lead(inp, gamma_dis = 1.6), 30)
  expect_equal(inc[1], 0.0029)
  expect_true(all(diff(inc) > 0))
  expect_equal(p_lead_dis(inp, 5), rep(0, 5))   # p_dis_yr = 0
})

test_that("t_lead_cap flattens the lead hazard after the cap", {
  x <- with_lead(inp, gamma_dis = 1.6)
  x$t_lead_cap <- 10
  p <- p_lead_dis(x, 30)
  uncapped <- p_lead_dis(with_lead(inp, gamma_dis = 1.6), 30)
  expect_equal(p[1:10], uncapped[1:10])
  expect_equal(p[11:30], rep(uncapped[10], 20))
})

test_that("late lead costs and disutilities apply by cycle in ICD and shock states only", {
  x <- with_lead(inp)
  J <- 30
  T <- x$t_repl
  p_dis <- p_lead_dis(x, J)

  cost0 <- rc_unit(J, inp, "icd")
  cost1 <- rc_unit(J, x, "icd")
  expected_extra <- x$p_inf_yr * x$c_inf + p_dis * x$c_dis
  # ICD tunnel state 2 and shock tunnel state, all cycles
  expect_equal(cost1[2, ] - cost0[2, ], expected_extra)
  expect_equal(cost1[T + 1, ] - cost0[T + 1, ], expected_extra)
  # death, no-ICD and SCD states untouched
  expect_equal(cost1[(2 * T + 1):(2 * T + 4), ], cost0[(2 * T + 1):(2 * T + 4), ])

  u0 <- re_unit(J, inp, "icd")
  u1 <- re_unit(J, x, "icd")
  expect_equal(u1[2, ] - u0[2, ], (x$p_inf_yr + p_dis) * x$u_compl)
  expect_equal(u1[(2 * T + 1):(2 * T + 4), ], u0[(2 * T + 1):(2 * T + 4), ])

  # no_icd strategy gets no late lead complications
  expect_equal(rc_unit(J, x, "no_icd"), rc_unit(J, inp, "no_icd"))
  expect_equal(re_unit(J, x, "no_icd"), re_unit(J, inp, "no_icd"))
})

test_that("lead age is not reset at generator replacement", {
  x <- with_lead(inp)
  J <- 30
  cost <- rc_unit(J, x, "icd") - rc_unit(J, inp, "icd")
  # state 1 at cycle 12 (post-replacement) carries the lead hazard of year 12 only,
  # plus the replacement complication cost
  extra_12 <- x$p_inf_yr * x$c_inf + p_lead_dis(x, J)[12] * x$c_dis
  expect_equal(cost[1, 12], extra_12)
})

test_that("replacement cost includes the replacement complication cost", {
  J <- 30
  mat <- rc_unit(J, inp, "icd")
  expect_equal(mat[1, 2], inp$c_appt + inp$c_icd_repl +
                 inp$p_repl_compl * c_compl_weighted(inp))
  expect_equal(mat[1, 1], inp$c_appt)
})

test_that("zero new probabilities leave old replacement behaviour except p_repl_compl", {
  x <- inp
  x$p_repl_compl <- 0
  mat <- rc_unit(30, x, "icd")
  expect_equal(mat[1, 2], x$c_appt + x$c_icd_repl)
  expect_equal(re_unit(30, x, "icd")[1, 2],
               x$q_hcm * x$u_icd + x$u_repl)
})

test_that("c_compl_weighted is zero when both initial probabilities are zero", {
  x <- inp
  x$p_inf_init <- 0
  x$p_dis_init <- 0
  expect_equal(c_compl_weighted(x), 0)
})

test_that("sample_param shares common random numbers across device labels", {
  dist <- readr::read_csv(
    system.file("extdata/input_distns.csv", package = "HCMSCDCEanalysis"),
    show_col_types = FALSE)
  tv <- dist[dist$label == "tvicd", ]
  sc <- dist[dist$label == "sicd", ]
  set.seed(42); a <- do.call(rbind, lapply(1:200, function(i) sample_param(tv)))
  set.seed(42); b <- do.call(rbind, lapply(1:200, function(i) sample_param(sc)))
  # parameters with identical distributions are identical draw-for-draw
  expect_equal(a$q_hcm, b$q_hcm)
  expect_equal(a$u_icd, b$u_icd)
  expect_equal(a$u_implant, b$u_implant)
  # the point and beta/normal rows follow each device's own distribution
  expect_true(all(a$t_repl == 10))
  expect_true(all(b$t_repl == 8))
  expect_gt(mean(a$p_compl), mean(b$p_compl))
})

test_that("qpoint returns the fixed value", {
  expect_equal(qpoint(c(0.1, 0.9), 3), c(3, 3))
})

test_that("beta probability rows sample within (0, 1) around their mean", {
  dist <- readr::read_csv(
    system.file("extdata/input_distns.csv", package = "HCMSCDCEanalysis"),
    show_col_types = FALSE)
  means <- list(
    p_compl      = c(tvicd = 0.043,  sicd = 0.015),
    p_repl_compl = c(tvicd = 0.0091, sicd = 0.002),
    p_inf_yr     = c(tvicd = 0.002,  sicd = 0.0005),
    p_dis_yr     = c(tvicd = 0.0029, sicd = 0.0027)
  )
  for (v in names(means)) {
    for (lbl in c("tvicd", "sicd")) {
      d <- dist[dist$label == lbl & dist$variable == v, ]
      expect_equal(d$distn, "beta")
      set.seed(1)
      draws <- replicate(2000, sample_param(d)[[v]])
      expect_true(all(draws >= 0 & draws < 1))
      expect_equal(mean(draws), means[[v]][[lbl]], tolerance = 0.15)
    }
  }
})

test_that("S-ICD has lower complication inputs than TV-ICD", {
  vals <- readr::read_csv(
    system.file("extdata/input_values.csv", package = "HCMSCDCEanalysis"),
    show_col_types = FALSE)
  get <- function(v, lbl) vals$value[vals$variable == v & vals$label == lbl]
  for (v in c("p_compl", "p_repl_compl", "p_inf_yr", "p_dis_yr")) {
    expect_lt(get(v, "sicd"), get(v, "tvicd"))
  }
})
