# Tests for input validation

# Skip all tests on CRAN to avoid timeouts
skip_on_cran()

# ==============================================================================
# Y VALIDATION TESTS
# ==============================================================================

test_that("lpmec_onerun errors on NULL Y", {
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))
  expect_error(lpmec_onerun(NULL, obs), "'Y' is required")
})

test_that("lpmec_onerun errors on non-numeric Y", {
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))
  Y_char <- rep("a", 80)
  expect_error(lpmec_onerun(Y_char, obs), "'Y' must be a numeric vector")
})

test_that("lpmec_onerun errors on Y with too few observations", {
  obs <- as.data.frame(matrix(sample(c(0, 1), 5 * 6, replace = TRUE), ncol = 6))
  Y <- rnorm(5)
  expect_error(lpmec_onerun(Y, obs), "at least 10 observations")
})

test_that("lpmec_onerun errors on all-NA Y", {
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))
  Y <- rep(NA_real_, 80)
  expect_error(lpmec_onerun(Y, obs), "cannot be all NA")
})

# ==============================================================================
# OBSERVABLES VALIDATION TESTS
# ==============================================================================

test_that("lpmec_onerun errors on NULL observables", {
  Y <- rnorm(80)
  expect_error(lpmec_onerun(Y, NULL), "'observables' is required")
})

test_that("lpmec_onerun errors on observables with wrong number of rows", {
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 50 * 6, replace = TRUE), ncol = 6))
  expect_error(lpmec_onerun(Y, obs), "must match length of 'Y'")
})

test_that("lpmec_onerun errors on observables with too few columns", {
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 3, replace = TRUE), ncol = 3))
  expect_error(lpmec_onerun(Y, obs), "at least 4 columns")
})

test_that("lpmec_onerun warns on excessive missing data", {
  set.seed(456)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 8, replace = TRUE), ncol = 8))

  # Set > 50% of values to NA in a controlled way
  # Set columns 5-8 completely to NA (50% of columns), plus some extra
  obs[, 5:8] <- NA
  obs[1:10, 1] <- NA  # Add extra NA to push over 50%

  # The warning should be raised before any computation starts
  # Use tryCatch to verify warning is raised even if function errors later
  warning_raised <- FALSE
  tryCatch(
    withCallingHandlers(
      lpmec_onerun(Y, obs, estimation_method = "averaging"),
      warning = function(w) {
        if (grepl("More than 50%", conditionMessage(w))) {
          warning_raised <<- TRUE
        }
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL  # Ignore errors, we only care about the warning
  )
  expect_true(warning_raised)
})

# ==============================================================================
# ESTIMATION METHOD VALIDATION TESTS
# ==============================================================================

test_that("lpmec_onerun errors on invalid estimation_method", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))
  expect_error(
    lpmec_onerun(Y, obs, estimation_method = "invalid_method"),
    "'estimation_method' must be one of"
  )
})

test_that("lpmec_onerun errors when custom method lacks function", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))
  expect_error(
    lpmec_onerun(Y, obs, estimation_method = "custom"),
    "'latent_estimation_fn' is required"
  )
})

test_that("lpmec_onerun errors when latent_estimation_fn is not a function", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))
  expect_error(
    lpmec_onerun(Y, obs, estimation_method = "custom", latent_estimation_fn = "not_a_function"),
    "'latent_estimation_fn' must be a function"
  )
})

# ==============================================================================
# ORDINAL PARAMETER VALIDATION TESTS
# ==============================================================================

test_that("lpmec_onerun errors on invalid ordinal parameter", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))
  expect_error(
    lpmec_onerun(Y, obs, estimation_method = "pca", ordinal = "TRUE"),
    "'ordinal' must be a single logical"
  )
})

# ==============================================================================
# MCMC_CONTROL VALIDATION TESTS
# ==============================================================================

test_that("lpmec_onerun errors on invalid mcmc_control backend", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))
  expect_error(
    lpmec_onerun(Y, obs, estimation_method = "pca",
                mcmc_control = list(backend = "invalid_backend")),
    "mcmc_control\\$backend must be either"
  )
})

test_that("lpmec_onerun errors on invalid mcmc_control n_samples_warmup", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))
  expect_error(
    lpmec_onerun(Y, obs, estimation_method = "pca",
                mcmc_control = list(n_samples_warmup = -1)),
    "n_samples_warmup must be a positive integer"
  )
})

# ==============================================================================
# LPME-SPECIFIC VALIDATION TESTS
# ==============================================================================

test_that("lpmec errors on invalid n_boot", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))
  expect_error(
    lpmec(Y, obs, n_boot = -1, n_partition = 1, estimation_method = "pca"),
    "'n_boot' must be a single positive integer"
  )
})

test_that("lpmec errors on invalid n_partition", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))
  expect_error(
    lpmec(Y, obs, n_boot = 1, n_partition = 0, estimation_method = "pca"),
    "'n_partition' must be a single positive integer"
  )
})

test_that("lpmec errors on boot_basis with wrong length", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))
  expect_error(
    lpmec(Y, obs, n_boot = 1, n_partition = 1, boot_basis = 1:50,
         estimation_method = "pca"),
    "'boot_basis' must have the same length"
  )
})

test_that("lpmec errors on invalid return_intermediaries", {
  set.seed(123)
  Y <- rnorm(80)
  obs <- as.data.frame(matrix(sample(c(0, 1), 80 * 6, replace = TRUE), ncol = 6))
  expect_error(
    lpmec(Y, obs, n_boot = 1, n_partition = 1, return_intermediaries = "yes",
         estimation_method = "pca"),
    "'return_intermediaries' must be a single logical"
  )
})
