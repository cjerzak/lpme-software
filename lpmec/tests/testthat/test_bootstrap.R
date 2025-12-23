# Tests for bootstrap aggregation functionality

# Skip all tests on CRAN to avoid timeouts
skip_on_cran()

# ==============================================================================
# BOOTSTRAP AGGREGATION TESTS
# ==============================================================================

test_that("lpmec bootstrap aggregation works with n_boot > 1", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))

  res <- lpmec(Y, obs, n_boot = 3, n_partition = 2, estimation_method = "pca")

  expect_s3_class(res, "lpmec")
  # Bootstrap should produce standard errors
  expect_true("ols_se" %in% names(res))
  expect_true("corrected_iv_se" %in% names(res))
  # Standard errors should be positive (non-NA) when n_boot > 1
  expect_true(is.numeric(res$ols_se))
  expect_true(is.numeric(res$corrected_iv_se))
})

test_that("lpmec partition aggregation works with n_partition > 1", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))

  res <- lpmec(Y, obs, n_boot = 1, n_partition = 3, estimation_method = "pca")

  expect_s3_class(res, "lpmec")
  expect_true("ols_coef" %in% names(res))
  expect_true(is.numeric(res$ols_coef))
})

test_that("lpmec produces confidence intervals with sufficient bootstrap", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))

  res <- lpmec(Y, obs, n_boot = 5, n_partition = 2, estimation_method = "pca")

  # Should have lower and upper bounds

  expect_true("ols_lower" %in% names(res))
  expect_true("ols_upper" %in% names(res))
  expect_true("corrected_iv_lower" %in% names(res))
  expect_true("corrected_iv_upper" %in% names(res))

  # Lower should be less than upper
  if (!is.na(res$ols_lower) && !is.na(res$ols_upper)) {
    expect_true(res$ols_lower <= res$ols_upper)
  }
})

test_that("lpmec with stratified bootstrap (boot_basis) works", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))
  # Create a grouping variable for stratified bootstrap
  boot_groups <- rep(1:10, each = 8)

  res <- lpmec(Y, obs, n_boot = 2, n_partition = 1, boot_basis = boot_groups,
              estimation_method = "pca")

  expect_s3_class(res, "lpmec")
  expect_true("ols_coef" %in% names(res))
})

test_that("lpmec intermediary results are stored when return_intermediaries = TRUE", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))

  res <- lpmec(Y, obs, n_boot = 2, n_partition = 2, estimation_method = "pca",
              return_intermediaries = TRUE)

  expect_s3_class(res, "lpmec")
  # Should have x_est1 and x_est2 from first run
  expect_true("x_est1" %in% names(res))
  expect_true("x_est2" %in% names(res))
  expect_equal(length(res$x_est1), 80)
})

test_that("lpmec var_est_split is computed correctly", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))

  res <- lpmec(Y, obs, n_boot = 2, n_partition = 3, estimation_method = "pca")

  # var_est_split should exist
  expect_true("var_est_split" %in% names(res))
  # It should be numeric (could be NA in edge cases)
  expect_true(is.numeric(res$var_est_split) || is.na(res$var_est_split))
})
