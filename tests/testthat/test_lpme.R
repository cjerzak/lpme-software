test_that("lpme basic functionality", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0,1), 80*4, replace=TRUE), ncol=4))
  res <- lpme(Y, obs, n_boot = 1, n_partition = 1, estimation_method = "pca")
  expect_s3_class(res, "lpme")
  expect_true("ols_coef" %in% names(res))
})

test_that("lpme orientation_signs validation", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0,1), 80*4, replace=TRUE), ncol=4))
  expect_error(
    lpme(Y, obs, n_boot = 1, n_partition = 1, estimation_method = "pca",
         orientation_signs = c(1, 0, 1)),
    "orientation_signs"
  )
})

test_that("lpme S3 methods", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0,1), 80*4, replace=TRUE), ncol=4))
  res <- lpme(Y, obs, n_boot = 1, n_partition = 1, estimation_method = "pca")
  sum_df <- summary(res)
  expect_s3_class(sum_df, "data.frame")
  expect_equal(
    row.names(sum_df),
    c("OLS", "IV", "Corrected IV", "Corrected OLS", "Bayesian OLS (Outer)", "Bayesian OLS (Inner)")
  )
  expect_type(capture.output(print(res)), "character")
  expect_silent(plot(res))
  expect_error(plot(res, type = "coefficients"))
})

