test_that("inits_by_interv works", {
  # Create a mock data set
  mock_data <- data.frame(
    icd = c(TRUE, FALSE, TRUE),
    risk_over_6 = c(FALSE, TRUE, FALSE),
    risk_over_4 = c(TRUE, TRUE, FALSE)
  )

  res <- inits_by_interv(mock_data, c("icd", "risk_over_6", "risk_over_4"))
  expect_s3_class(res, "data.frame")
  expect_true("risk" %in% names(res))
  expect_true("status" %in% names(res))
  expect_true("n" %in% names(res))
})
