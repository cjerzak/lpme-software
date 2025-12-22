# Tests using the included KnowledgeVoteDuty dataset

# Skip all tests on CRAN to avoid timeouts
skip_on_cran()

# ==============================================================================
# KNOWLEDGEVOTEDUTY DATASET TESTS
# ==============================================================================

test_that("lpme works with KnowledgeVoteDuty dataset using PCA", {
  skip_if_not_installed("lpme")

  data(KnowledgeVoteDuty, package = "lpme")

  # Extract outcome and observables
  Y <- KnowledgeVoteDuty$voteduty
  obs <- KnowledgeVoteDuty[, c("SenateTerm", "SpendLeast", "HouseParty", "SenateParty")]

  # Run with minimal bootstrap for speed
  res <- lpme(Y, obs, n_boot = 1, n_partition = 1, estimation_method = "pca")

  expect_s3_class(res, "lpme")
  expect_true("ols_coef" %in% names(res))
  expect_true("corrected_iv_coef" %in% names(res))
  expect_true(is.numeric(res$ols_coef))
  expect_true(is.numeric(res$corrected_iv_coef))
})

test_that("lpme_onerun works with KnowledgeVoteDuty dataset", {
  skip_if_not_installed("lpme")

  data(KnowledgeVoteDuty, package = "lpme")

  Y <- KnowledgeVoteDuty$voteduty
  obs <- KnowledgeVoteDuty[, c("SenateTerm", "SpendLeast", "HouseParty", "SenateParty")]

  res <- lpme_onerun(Y, obs, estimation_method = "pca")

  expect_s3_class(res, "lpme_onerun")
  expect_equal(length(res$x_est1), nrow(KnowledgeVoteDuty))
  expect_equal(length(res$x_est2), nrow(KnowledgeVoteDuty))
})

test_that("lpme with KnowledgeVoteDuty produces sensible coefficient directions", {
  skip_if_not_installed("lpme")

  data(KnowledgeVoteDuty, package = "lpme")

  Y <- KnowledgeVoteDuty$voteduty
  obs <- KnowledgeVoteDuty[, c("SenateTerm", "SpendLeast", "HouseParty", "SenateParty")]

  res <- lpme_onerun(Y, obs, estimation_method = "pca")

  # Latent variables should be correlated with row means of observables
  knowledge_simple <- rowMeans(obs)
  cor_est1 <- cor(res$x_est1, knowledge_simple, use = "complete.obs")
  cor_est2 <- cor(res$x_est2, knowledge_simple, use = "complete.obs")

  # Should have positive correlation with simple average (same construct)
  expect_true(abs(cor_est1) > 0.5)
  expect_true(abs(cor_est2) > 0.5)
})

test_that("lpme with KnowledgeVoteDuty using EM estimation", {
  skip_if_not_installed("lpme")
  skip_if_not_installed("emIRT")

  data(KnowledgeVoteDuty, package = "lpme")

  Y <- KnowledgeVoteDuty$voteduty
  obs <- KnowledgeVoteDuty[, c("SenateTerm", "SpendLeast", "HouseParty", "SenateParty")]

  res <- lpme(Y, obs, n_boot = 1, n_partition = 1, estimation_method = "em")

  expect_s3_class(res, "lpme")
  expect_true(is.numeric(res$ols_coef))
})

test_that("lpme with KnowledgeVoteDuty using averaging estimation", {
  skip_if_not_installed("lpme")

  data(KnowledgeVoteDuty, package = "lpme")

  Y <- KnowledgeVoteDuty$voteduty
  obs <- KnowledgeVoteDuty[, c("SenateTerm", "SpendLeast", "HouseParty", "SenateParty")]

  res <- lpme(Y, obs, n_boot = 1, n_partition = 1, estimation_method = "averaging")

  expect_s3_class(res, "lpme")
  expect_true(is.numeric(res$ols_coef))
})

test_that("KnowledgeVoteDuty dataset has expected structure", {
  skip_if_not_installed("lpme")

  data(KnowledgeVoteDuty, package = "lpme")

  # Check it's a data frame

  expect_s3_class(KnowledgeVoteDuty, "data.frame")

  # Check expected columns exist
  expected_cols <- c("voteduty", "SenateTerm", "SpendLeast", "HouseParty", "SenateParty")
  expect_true(all(expected_cols %in% colnames(KnowledgeVoteDuty)))

  # Check dimensions (should have ~3059 rows per documentation)
  expect_true(nrow(KnowledgeVoteDuty) > 3000)
  expect_equal(ncol(KnowledgeVoteDuty), 5)

  # Check binary columns are actually binary
  for (col in c("SenateTerm", "SpendLeast", "HouseParty", "SenateParty")) {
    expect_true(all(KnowledgeVoteDuty[[col]] %in% c(0, 1, NA)))
  }

  # Check voteduty range (1-7 per documentation)
  expect_true(all(KnowledgeVoteDuty$voteduty >= 1 & KnowledgeVoteDuty$voteduty <= 7,
                  na.rm = TRUE))
})
