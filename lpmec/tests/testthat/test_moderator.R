# Skip all tests on CRAN to avoid timeouts
skip_on_cran()

# Small helper: a moderator experiment with M pre-computed noisy scores of the
# same latent X (scores mode; no item batteries).
make_scores_moderator_data <- function(n = 600L,
                                       n_scores = 3L,
                                       score_sd = 0.7,
                                       b_treatment = 0.2,
                                       b_moderator = 0.2,
                                       b_interaction = 0.3,
                                       seed = 404L) {
  set.seed(seed)
  X <- rnorm(n)
  treatment <- rbinom(n, 1, 0.5)
  Y <- b_treatment * treatment + b_moderator * X +
    b_interaction * treatment * X + rnorm(n)
  scores <- lapply(seq_len(n_scores), function(j) X + rnorm(n, sd = score_sd))
  names(scores) <- paste0("s", seq_len(n_scores))
  list(Y = Y, treatment = treatment, X = X, scores = scores)
}

test_that(".lpmec_fit_interaction matches lm(Y ~ T*x) to 1e-10", {
  set.seed(11)
  n <- 200
  tr <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 0.2 * tr + 0.3 * x + 0.25 * tr * x + rnorm(n)

  fit <- .lpmec_fit_interaction(y, tr, x)
  lm_fit <- stats::lm(y ~ tr * x)
  lm_coef <- stats::coef(lm_fit)
  lm_se <- summary(lm_fit)$coefficients[, "Std. Error"]

  expect_equal(fit$treatment_coef, unname(lm_coef[["tr"]]), tolerance = 1e-10)
  expect_equal(fit$main_coef, unname(lm_coef[["x"]]), tolerance = 1e-10)
  expect_equal(fit$interaction_coef, unname(lm_coef[["tr:x"]]),
               tolerance = 1e-10)
  expect_equal(fit$treatment_se, unname(lm_se[["tr"]]), tolerance = 1e-10)
  expect_equal(fit$main_se, unname(lm_se[["x"]]), tolerance = 1e-10)
  expect_equal(fit$interaction_se, unname(lm_se[["tr:x"]]), tolerance = 1e-10)
  expect_equal(fit$n_obs, n)

  # with covariates
  covariate_matrix <- cbind(c1 = rnorm(n), c2 = runif(n))
  y2 <- y + 0.4 * covariate_matrix[, "c1"] - 0.2 * covariate_matrix[, "c2"]
  fit2 <- .lpmec_fit_interaction(y2, tr, x, covariate_matrix)
  lm_fit2 <- stats::lm(
    y2 ~ tr * x + c1 + c2,
    data = data.frame(y2 = y2, tr = tr, x = x, covariate_matrix)
  )
  lm_coef2 <- stats::coef(lm_fit2)
  lm_se2 <- summary(lm_fit2)$coefficients[, "Std. Error"]

  expect_equal(fit2$interaction_coef, unname(lm_coef2[["tr:x"]]),
               tolerance = 1e-10)
  expect_equal(fit2$interaction_se, unname(lm_se2[["tr:x"]]),
               tolerance = 1e-10)
  expect_equal(fit2$main_coef, unname(lm_coef2[["x"]]), tolerance = 1e-10)
  expect_equal(unname(fit2$coef_all[c("c1", "c2")]),
               unname(lm_coef2[c("c1", "c2")]), tolerance = 1e-10)
  expect_equal(unname(fit2$se_all[c("c1", "c2")]),
               unname(lm_se2[c("c1", "c2")]), tolerance = 1e-10)

  # rows with non-finite entries are dropped (as lm drops NAs)
  x_na <- x
  x_na[1:5] <- NA_real_
  fit3 <- .lpmec_fit_interaction(y, tr, x_na)
  lm_fit3 <- stats::lm(y ~ tr * x_na)
  expect_equal(fit3$interaction_coef,
               unname(stats::coef(lm_fit3)[["tr:x_na"]]), tolerance = 1e-10)
  expect_equal(fit3$n_obs, n - 5L)
})

test_that("Spearman-Brown step-up is exact in items mode (M = 1)", {
  dat <- make_moderator_test_data(n = 400L, seed = 205L)
  set.seed(1)
  run <- lpmec_moderator_onerun(
    Y = dat$Y,
    treatment = dat$treatment,
    observables = list(knowledge = dat$items),
    estimation_method = "averaging"
  )
  expect_s3_class(run, "lpmec_moderator_onerun")
  expect_equal(run$n_measures, 1L)
  expect_equal(run$sb_factor, 2)

  r_half <- stats::cor(run$x_est1[, 1L], run$x_est2[, 1L],
                       use = "pairwise.complete.obs")
  expect_equal(run$rho_half, r_half, tolerance = 1e-12)
  expect_equal(run$rho_score, 2 * r_half / (1 + r_half), tolerance = 1e-12)
  expect_equal(run$correction_factor, sqrt(run$rho_score), tolerance = 1e-12)
  expect_equal(run$corrected_interaction_coef,
               run$interaction_coef / sqrt(run$rho_score), tolerance = 1e-12)
  expect_equal(run$corrected_main_coef,
               run$main_coef / sqrt(run$rho_score), tolerance = 1e-12)

  # the regressor is the z-scored full-battery score
  expect_equal(run$x_used, unname(run$x_est[, 1L]), tolerance = 1e-12)
  expect_null(run$per_measure)
  expect_null(run$cor_matrix)
})

test_that("generalized Spearman-Brown is exact in scores mode (M = 3)", {
  dat <- make_scores_moderator_data(n = 600L, n_scores = 3L, seed = 405L)
  run <- lpmec_moderator_onerun(
    Y = dat$Y,
    treatment = dat$treatment,
    scores = dat$scores
  )
  expect_equal(run$n_measures, 3L)
  expect_equal(run$sb_factor, 3)

  cor_matrix <- stats::cor(run$x_est, use = "pairwise.complete.obs")
  r_bar <- mean(cor_matrix[upper.tri(cor_matrix)])
  expect_equal(run$rho_half, r_bar, tolerance = 1e-12)
  expect_equal(run$rho_score, 3 * r_bar / (1 + 2 * r_bar), tolerance = 1e-12)
  expect_equal(run$corrected_interaction_coef,
               run$interaction_coef / sqrt(run$rho_score), tolerance = 1e-12)

  # regressor is the z-score of the mean of the aligned z-scored measures
  row_mean <- rowMeans(run$x_est)
  expect_equal(run$x_used,
               (row_mean - mean(row_mean)) / stats::sd(row_mean),
               tolerance = 1e-12)
})

test_that("corrected interaction recovers the plim at n = 4000", {
  dat <- make_moderator_test_data(n = 4000L, seed = 202L)
  set.seed(2)
  run <- lpmec_moderator_onerun(
    Y = dat$Y,
    treatment = dat$treatment,
    observables = list(knowledge = dat$items),
    estimation_method = "averaging"
  )

  # naive is attenuated by sqrt(rho); corrected recovers b_interaction
  expect_lt(abs(run$corrected_interaction_coef - dat$b_interaction), 0.08)
  expect_lt(abs(run$interaction_coef -
                  dat$b_interaction * sqrt(run$rho_score)), 0.08)
  expect_lt(run$interaction_coef, run$corrected_interaction_coef)
  expect_gt(run$rho_score, run$rho_half)  # step-up raises the half reliability
  expect_lt(run$rho_score, 1)

  # the main-effect correction follows the same algebra
  expect_equal(run$corrected_main_coef,
               run$main_coef / sqrt(run$rho_score), tolerance = 1e-12)
})

test_that("whole-pipeline bootstrap returns the expected shapes", {
  dat <- make_moderator_test_data(n = 250L, seed = 303L)
  result <- suppressMessages(lpmec_moderator(
    Y = dat$Y,
    treatment = dat$treatment,
    observables = list(knowledge = dat$items),
    estimation_method = "averaging",
    n_boot = 2L,
    n_partition = 2L,
    seed = 9
  ))
  expect_s3_class(result, "lpmec_moderator")

  # (n_boot + 1) x n_partition runs, in order
  expect_equal(result$Intermediary_BootIndex, rep(1:3, each = 2L))
  expect_equal(result$Intermediary_PartitionIndex, rep(1:2, times = 3L))
  for (field in c("interaction_coef", "corrected_interaction_coef",
                  "rho_half", "rho_score")) {
    expect_length(result[[paste0("Intermediary_", field)]], 6L)
  }

  # headline quantities with bootstrap summaries, including the reliabilities
  for (field in c("interaction", "corrected_interaction", "rho_half",
                  "rho_score")) {
    prefix <- if (grepl("^rho", field)) field else paste0(field, "_coef")
    expect_true(is.finite(result[[prefix]]))
    expect_true(is.finite(result[[paste0(field, "_se")]]))
    expect_lte(result[[paste0(field, "_lower")]],
               result[[paste0(field, "_upper")]])
  }
  expect_equal(result$n_boot, 2L)
  expect_equal(result$n_partition, 2L)
  expect_equal(result$n_obs, 250L)

  # whole pipeline is re-run per (boot, partition): fresh splits give fresh
  # reliability estimates across runs
  expect_gt(length(unique(result$Intermediary_rho_half)), 1L)

  # reproducible under the same seed
  result2 <- suppressMessages(lpmec_moderator(
    Y = dat$Y,
    treatment = dat$treatment,
    observables = list(knowledge = dat$items),
    estimation_method = "averaging",
    n_boot = 2L,
    n_partition = 2L,
    seed = 9
  ))
  expect_equal(result2$corrected_interaction_coef,
               result$corrected_interaction_coef)
  expect_equal(result2$Intermediary_rho_score, result$Intermediary_rho_score)

  # stratified bootstrap on the treatment arms runs and keeps shapes
  strat <- suppressMessages(lpmec_moderator(
    Y = dat$Y,
    treatment = dat$treatment,
    observables = list(knowledge = dat$items),
    estimation_method = "averaging",
    n_boot = 1L,
    n_partition = 1L,
    boot_basis = dat$treatment,
    seed = 10
  ))
  expect_equal(strat$Intermediary_BootIndex, 1:2)
  expect_true(is.finite(strat$corrected_interaction_coef))
})

test_that("n_partition is coerced to 1 with a message outside items mode", {
  dat <- make_scores_moderator_data(n = 300L, n_scores = 2L, seed = 406L)
  messages_seen <- testthat::capture_messages(
    result <- lpmec_moderator(
      Y = dat$Y,
      treatment = dat$treatment,
      scores = dat$scores,
      n_boot = 0L,
      n_partition = 3L
    )
  )
  expect_true(any(grepl("n_partition", messages_seen)))
  expect_equal(result$n_partition, 1L)
  expect_equal(result$Intermediary_BootIndex, 1L)

  # without bootstrap, uncertainty summaries are NA
  expect_true(is.na(result$corrected_interaction_se))
  expect_true(is.na(result$rho_score_lower))
})

test_that("reliability floors yield NA corrections with one warning", {
  set.seed(77)
  n <- 300
  treatment <- rbinom(n, 1, 0.5)
  X <- rnorm(n)
  Y <- 0.2 * treatment + 0.3 * X + rnorm(n)

  # two exactly orthogonal noise scores: rho_half = 0 (sign alignment cannot
  # rescue a zero correlation), so rho_score = 0 < min_reliability
  noise_a <- rnorm(n)
  noise_b <- stats::residuals(stats::lm(rnorm(n) ~ noise_a))
  scores_bad <- list(a = noise_a, b = as.numeric(noise_b))
  warnings_seen <- testthat::capture_warnings(
    run <- lpmec_moderator_onerun(Y, treatment, scores = scores_bad)
  )
  expect_length(warnings_seen, 1L)
  expect_match(warnings_seen, "min_reliability")
  expect_true(run$rho_score < 0.05 || is.na(run$rho_score))
  expect_true(is.na(run$corrected_interaction_coef))
  expect_true(is.na(run$corrected_main_coef))
  expect_true(is.na(run$correction_factor))
  expect_true(run$reliability_floored)
  expect_true(is.finite(run$interaction_coef))  # naive is still reported

  # a single scores-mode measure has no reliability information: rho_score is
  # non-finite -> NA + one warning
  warnings_single <- testthat::capture_warnings(
    run_single <- lpmec_moderator_onerun(
      Y, treatment, scores = list(a = X + rnorm(n))
    )
  )
  expect_length(warnings_single, 1L)
  expect_match(warnings_single, "min_reliability")
  expect_true(is.na(run_single$rho_half))
  expect_true(is.na(run_single$rho_score))
  expect_true(is.na(run_single$corrected_interaction_coef))
})

test_that("treatment validation errors and continuous-treatment message", {
  dat <- make_scores_moderator_data(n = 200L, n_scores = 2L, seed = 407L)
  Y <- dat$Y
  scores <- dat$scores
  n <- length(Y)

  expect_error(
    lpmec_moderator_onerun(Y, scores = scores),
    "'treatment' is required"
  )
  expect_error(
    lpmec_moderator_onerun(Y, treatment = as.character(dat$treatment),
                           scores = scores),
    "numeric"
  )
  expect_error(
    lpmec_moderator_onerun(Y, treatment = dat$treatment[-1L], scores = scores),
    "length"
  )
  treatment_na <- dat$treatment
  treatment_na[3L] <- NA_real_
  expect_error(
    lpmec_moderator_onerun(Y, treatment = treatment_na, scores = scores),
    "finite"
  )
  expect_error(
    lpmec_moderator_onerun(Y, treatment = rep(1, n), scores = scores),
    "variance"
  )

  # more than 2 unique values: message, not error
  set.seed(88)
  dose <- rnorm(n)
  expect_message(
    run <- lpmec_moderator_onerun(Y, treatment = dose, scores = scores),
    "unique"
  )
  expect_true(is.finite(run$interaction_coef))

  # the aggregator validates the treatment up front too
  expect_error(
    lpmec_moderator(Y, treatment = rep(0, n), scores = scores, n_boot = 0L),
    "variance"
  )
})

test_that("per_measure is present iff M >= 3", {
  dat3 <- make_scores_moderator_data(n = 500L, n_scores = 3L, seed = 408L)
  run3 <- lpmec_moderator_onerun(
    Y = dat3$Y, treatment = dat3$treatment, scores = dat3$scores
  )
  expect_s3_class(run3$per_measure, "data.frame")
  expect_equal(nrow(run3$per_measure), 3L)
  expect_equal(run3$per_measure$measure, c("s1", "s2", "s3"))
  expect_equal(names(run3$per_measure),
               c("measure", "interaction_coef", "interaction_se", "triad",
                 "corrected_interaction_coef"))
  expect_true(all(is.finite(run3$per_measure$triad)))
  expect_equal(run3$per_measure$corrected_interaction_coef,
               run3$per_measure$interaction_coef /
                 sqrt(run3$per_measure$triad),
               tolerance = 1e-12)
  # triads match the pooled cross-measure correlation matrix identity
  expect_equal(unname(run3$per_measure$triad),
               unname(.lpmec_triad_reliabilities(run3$cor_matrix)),
               tolerance = 1e-12)
  expect_equal(dim(run3$cor_matrix), c(3L, 3L))

  dat2 <- make_scores_moderator_data(n = 500L, n_scores = 2L, seed = 409L)
  run2 <- lpmec_moderator_onerun(
    Y = dat2$Y, treatment = dat2$treatment, scores = dat2$scores
  )
  expect_null(run2$per_measure)
  expect_equal(run2$sb_factor, 2)

  # M = 1 items mode also has no per-measure table (covered above);
  # M >= 3 survives the aggregator: point run's table is carried through
  agg3 <- suppressMessages(lpmec_moderator(
    Y = dat3$Y, treatment = dat3$treatment, scores = dat3$scores,
    n_boot = 1L, n_partition = 1L, seed = 12
  ))
  expect_s3_class(agg3$per_measure, "data.frame")
  expect_equal(nrow(agg3$per_measure), 3L)
})

test_that("covariates flow through the moderator pipeline", {
  dat <- make_scores_moderator_data(n = 400L, n_scores = 2L, seed = 410L)
  set.seed(13)
  covariate_df <- data.frame(age = rnorm(400L), female = rbinom(400L, 1, 0.5))
  run <- lpmec_moderator_onerun(
    Y = dat$Y,
    treatment = dat$treatment,
    scores = dat$scores,
    covariates = covariate_df
  )
  expect_equal(run$covariate_names, c("age", "female"))
  expect_true(all(c("age", "female") %in% names(run$coef_all)))
  expect_equal(run$corrected_interaction_coef,
               run$interaction_coef / sqrt(run$rho_score), tolerance = 1e-12)
})
