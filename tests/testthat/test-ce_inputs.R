test_that("lambda_tunnel works and builds matrices correctly", {
  # Create a mock 6x6 transition matrix
  # States: 1: ICD, 2: Shock, 3: Death ICD, 4: HCM (NoICD), 5: SCD, 6: Death NoICD
  lambda <- matrix(0, 6, 6)
  lambda[1, 1] <- 0.8  # ICD -> ICD
  lambda[1, 2] <- 0.1  # ICD -> Shock
  lambda[1, 3] <- 0.1  # ICD -> Death ICD
  lambda[4, 4] <- 0.9  # HCM -> HCM
  lambda[4, 5] <- 0.05 # HCM -> SCD
  lambda[4, 6] <- 0.05 # HCM -> Death NoICD
  
  # 1. Test T = 5 (n_tunnel = 5)
  n_tunnel <- 5
  res_mat <- lambda_tunnel(lambda, n_tunnel = n_tunnel, alpha = 0.2)
  
  # States should be 2 * 5 + 4 = 14
  expect_equal(dim(res_mat), c(14, 14))
  
  # All rows should sum to 1 (except maybe dead states which are sink states with diagonal 1)
  expect_equal(rowSums(res_mat), rep(1, 14))
  
  # Check ICD tunnel transitions (1 to 5):
  # ICD 1 (state 1) -> ICD 2 (state 2) with prob 0.8
  expect_equal(res_mat[1, 2], 0.8)
  # ICD 1 (state 1) -> Shock 1 (state 6) with prob 0.1
  expect_equal(res_mat[1, 6], 0.1)
  # ICD 1 (state 1) -> Death ICD (state 11) with prob 0.1
  expect_equal(res_mat[1, 11], 0.1)
  
  # ICD 5 (state 5) -> ICD 1 (state 1) (replaced) with prob 0.8
  expect_equal(res_mat[5, 1], 0.8)
  # ICD 5 (state 5) -> Shock 5 (state 10) with prob 0.1
  expect_equal(res_mat[5, 10], 0.1)
  
  # Check Shock tunnel transitions (6 to 10):
  # Shock 1 (state 6) -> ICD 2 (state 2) (progression, not skipping a year!) with prob 1 - pdeath = 0.9
  expect_equal(res_mat[6, 2], 0.9)
  expect_equal(res_mat[6, 11], 0.1) # Death
  
  # Shock 4 (state 9) -> ICD 5 (state 5) with prob 0.9
  expect_equal(res_mat[9, 5], 0.9)
  
  # Shock 5 (state 10) -> ICD 1 (state 1) (replaced) with prob 0.9
  expect_equal(res_mat[10, 1], 0.9)
  
  # Check HCM (state 12) transitions with alpha = 0.2
  # HCM -> ICD 1: alpha * lambda[4,4] = 0.2 * 0.9 = 0.18
  expect_equal(res_mat[12, 1], 0.18)
  # HCM -> HCM (state 12): (1 - alpha) * lambda[4,4] = 0.8 * 0.9 = 0.72
  expect_equal(res_mat[12, 12], 0.72)
  # HCM -> SCD (state 13): 0.05
  expect_equal(res_mat[12, 13], 0.05)
  # HCM -> Death NoICD (state 14): 0.05
  expect_equal(res_mat[12, 14], 0.05)

  # 2. Test T = 1 (n_tunnel = 1)
  res_mat_1 <- lambda_tunnel(lambda, n_tunnel = 1, alpha = 0)
  expect_equal(dim(res_mat_1), c(6, 6))
  expect_equal(rowSums(res_mat_1), rep(1, 6))
  # Shock 1 (state 2) -> ICD 1 (state 1) (replaced) with prob 1 - pdeath = 0.9
  expect_equal(res_mat_1[2, 1], 0.9)
})

test_that("get_init_counts works on ipd_risk data", {
  data("ipd_risk", package = "HCMSCDCEanalysis")
  
  counts <- get_init_counts(ipd_risk)
  expect_type(counts, "list")
  expect_setequal(names(counts), c("obs", "risk6", "risk4", "no_icd"))
  
  for (nm in names(counts)) {
    expect_true(is.numeric(counts[[nm]]$ICD))
    expect_true(is.numeric(counts[[nm]]$no_ICD))
  }
  
  # no_icd should have 0 ICD patients
  expect_equal(counts$no_icd$ICD, 0)
  expect_equal(counts$no_icd$no_ICD, nrow(ipd_risk))
})

test_that("create_pop_array creates correct populations", {
  init_counts <- list(
    obs = list(ICD = 10, no_ICD = 90),
    risk6 = list(ICD = 15, no_ICD = 85),
    risk4 = list(ICD = 20, no_ICD = 80),
    no_icd = list(ICD = 0, no_ICD = 100)
  )
  
  # Special case n_tunnel = 1 (6 states)
  pop_1 <- create_pop_array(init_counts, J = 10, n_sim = 2, n_tunnel = 1)
  expect_type(pop_1, "list")
  expect_setequal(names(pop_1), c("obs", "risk6", "risk4", "no_icd"))
  expect_equal(dim(pop_1$obs), c(6, 10, 2))
  expect_equal(pop_1$obs[, 1, 1], c(10, 0, 0, 90, 0, 0))
  
  # Tunnel model (n_tunnel = 5)
  pop_tunnel <- create_pop_array(init_counts, J = 10, n_sim = 2, n_tunnel = 5)
  expect_type(pop_tunnel, "list")
  expect_setequal(names(pop_tunnel), c("obs", "risk6", "risk4", "no_icd"))
  expect_equal(dim(pop_tunnel$obs), c(14, 10, 2))
  # State 1 is ICD 1. State 12 is NoICD Healthy (2 * 5 + 2 = 12)
  expect_equal(pop_tunnel$obs[1, 1, 1], 10)
  expect_equal(pop_tunnel$obs[12, 1, 1], 90)
  expect_equal(sum(pop_tunnel$obs[, 1, 1]), 100)
})

test_that("create_probs_array creates transition probability matrices", {
  # Mock list of lambda samples from MCMC (n_sim x 6 x 6)
  n_sim <- 3
  lambda_obs <- array(0, dim = c(n_sim, 6, 6))
  for (i in 1:n_sim) {
    lambda_obs[i, 1, 1] <- 0.8
    lambda_obs[i, 1, 2] <- 0.1
    lambda_obs[i, 1, 3] <- 0.1
    lambda_obs[i, 4, 4] <- 0.9
    lambda_obs[i, 4, 5] <- 0.05
    lambda_obs[i, 4, 6] <- 0.05
  }
  
  lambda_list <- list(
    obs = lambda_obs,
    risk6 = lambda_obs,
    risk4 = lambda_obs,
    no_icd = lambda_obs
  )
  
  # Special case n_tunnel = 1
  probs_1 <- create_probs_array(lambda_list, n_tunnel = 1)
  expect_equal(dim(probs_1$obs), c(6, 6, 1, n_sim))
  expect_equal(probs_1$obs[1, 2, 1, 1], 0.1) # ICD -> Shock
  
  # Tunnel model (n_tunnel = 5)
  probs_tunnel <- create_probs_array(lambda_list, n_tunnel = 5, alpha = 0.2)
  expect_equal(dim(probs_tunnel$obs), c(14, 14, 1, n_sim))
  expect_equal(probs_tunnel$obs[1, 2, 1, 1], 0.8) # ICD 1 -> ICD 2
  expect_equal(probs_tunnel$obs[6, 2, 1, 1], 0.9) # Shock 1 -> ICD 2
})

test_that("rc_unit applies replacement cost only in cycles j > 1 when n_tunnel > 1", {
  inp <- list(t_repl = 8, c_appt = 145, c_icd_repl = 17434, c_shock = 1766)
  mat <- rc_unit(J = 20, inp = inp, name_interv = "obs")
  expect_equal(dim(mat), c(20, 20))
  # Year 1 (j = 1) in State 1 should only have c_appt (initial implant is in rc_init)
  expect_equal(mat[1, 1], 145)
  # Years j > 1 in State 1 should have c_appt + c_icd_repl
  expect_equal(mat[1, 2:20], rep(145 + 17434, 19))
  # Other tunnel states (2..8) should only have c_appt across all years
  expect_equal(mat[2, ], rep(145, 20))
})
