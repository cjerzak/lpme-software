# Skip all tests on CRAN to avoid timeouts
skip_on_cran()

# ---------------------------------------------------------------------------
# Theory-consistency gates added after the ms/appendix audit:
#   1. split_scores mode: the split-based corrected OLS must use score-scale
#      (Spearman-Brown) reliabilities, recover beta, and agree with the
#      corrected split-IV (Prop 3 asymptotic OLS-IV agreement). The raw
#      half-split correction overcorrects and is rejected.
#   2. cross-measure IV at M >= 3: target-oriented triad correction
#      (Prop 3c), never sqrt of the pairwise correlation.
#   3. latent outcomes (Prop 4): every corrected estimator gains 1/sqrt(rho_Y).
#   4. design-local scale (Prop 2a): exact local-pooled bridge identity and
#      plim recovery of beta * sd(design-transformed X).
# ---------------------------------------------------------------------------

sb_up <- function(r) 2 * r / (1 + r)

test_that("split-scores corrected OLS recovers beta and agrees with split-IV", {
  dat <- make_panel_test_data(n_units = 500L, n_periods = 10L,
                              s = 0.5, sig2U = 0.5, beta = 0.4,
                              phi = 0.7, r2_within = 0.2, seed = 314L)
  run <- lpmec_panel_onerun(
    Y = dat$Y, unit = dat$unit, time = dat$time,
    split_scores = list(m = cbind(dat$t1, dat$t2)),
    design = "within"
  )

  # reported reliabilities: raw half-split correlations plus their step-up
  expect_equal(unname(run$split_rho_score[["m"]]),
               sb_up(run$split_correlation[["m"]]), tolerance = 1e-12)
  expect_equal(unname(run$design_split_rho_score[["m"]]),
               sb_up(run$design_split_correlation[["m"]]), tolerance = 1e-12)

  # correction algebra ties out against the stepped-up reliabilities
  expect_equal(
    run$corrected_ols_coef_split[["m"]],
    run$ols_coef[["m"]] * sqrt(run$split_rho_score[["m"]]) /
      run$design_split_rho_score[["m"]],
    tolerance = 1e-12
  )

  # plim recovery: the regressed score is the average of the two halves, so
  # the score-scale correction is consistent for beta ...
  expect_lt(abs(run$corrected_ols_coef[["m"]] - dat$beta), 0.05)
  expect_lt(abs(run$corrected_iv_coef[["m"]] - dat$beta), 0.05)
  # ... and corrected OLS agrees with corrected split-IV (Prop 3)
  expect_lt(abs(run$corrected_ols_coef[["m"]] - run$corrected_iv_coef[["m"]]),
            0.04)

  # the raw half-split correction (the pre-1.3.0 behavior) overcorrects
  raw_corrected <- run$ols_coef[["m"]] *
    sqrt(run$split_correlation[["m"]]) / run$design_split_correlation[["m"]]
  expect_gt(raw_corrected - dat$beta, 0.08)
})

test_that("cross-measure IV uses the target triad at M >= 3 (Prop 3c)", {
  dat <- make_panel_test_data(n_units = 500L, n_periods = 10L,
                              s = 0.5, sig2U = 0.5, beta = 0.4,
                              phi = 0.7, r2_within = 0.2, seed = 271L)
  set.seed(272L)
  n <- length(dat$X)
  s1 <- dat$X + rnorm(n, sd = sqrt(0.10))  # rho_1 = 0.909
  s2 <- dat$X + rnorm(n, sd = sqrt(0.80))  # rho_2 = 0.556
  s3 <- dat$X + rnorm(n, sd = sqrt(0.30))  # rho_3 = 0.769

  run <- lpmec_panel_onerun(
    Y = dat$Y, unit = dat$unit, time = dat$time,
    scores = list(m1 = s1, m2 = s2, m3 = s3),
    design = "within"
  )

  # structural identity: each oriented coefficient is multiplied by sqrt of
  # the REGRESSOR (target) measure's pooled triad reliability
  for (target in c("m1", "m2", "m3")) {
    for (instrument in setdiff(c("m1", "m2", "m3"), target)) {
      key <- paste0(target, "_by_", instrument)
      expect_equal(
        run$corrected_cross_iv_coef[[key]],
        run$cross_iv_coef[[key]] * sqrt(run$triad[[target]]),
        tolerance = 1e-12
      )
    }
  }

  # plim recovery in both orientations of the most unequal pair
  expect_lt(abs(run$corrected_cross_iv_coef[["m1_by_m2"]] - dat$beta), 0.06)
  expect_lt(abs(run$corrected_cross_iv_coef[["m2_by_m1"]] - dat$beta), 0.06)

  # the pairwise sqrt(r_ml) substitution is biased when reliabilities differ
  # (plim beta * (rho_l / rho_m)^(1/4)); it must not match the triad value
  pairwise_value <- run$cross_iv_coef[["m1_by_m2"]] *
    sqrt(run$cor_matrices[["pooled"]]["m1", "m2"])
  expect_gt(abs(pairwise_value - run$corrected_cross_iv_coef[["m1_by_m2"]]),
            0.02)
})

test_that("latent outcomes (Prop 4): pooled cross-section recovers beta", {
  set.seed(414L)
  n <- 4000L
  beta <- 0.5
  X <- rnorm(n)
  Y_latent <- beta * X + rnorm(n, sd = sqrt(1 - beta^2))  # Var(Y) = 1
  yh1 <- Y_latent + rnorm(n, sd = sqrt(0.4))
  yh2 <- Y_latent + rnorm(n, sd = sqrt(0.4))
  xh1 <- X + rnorm(n, sd = sqrt(0.5))
  xh2 <- X + rnorm(n, sd = sqrt(0.5))

  # unit may be NULL for the pooled design (pure cross-section)
  run <- lpmec_panel_onerun(
    Y_split_scores = cbind(yh1, yh2),
    split_scores = list(m = cbind(xh1, xh2)),
    design = "pooled"
  )

  # rho_Y: half correlation 1/1.4, stepped up to 1/1.2
  expect_equal(run$rho_Y_half, 1 / 1.4, tolerance = 0.03)
  expect_equal(run$rho_Y, sb_up(run$rho_Y_half), tolerance = 1e-12)

  # naive attenuated by sqrt(rho_X * rho_Y); corrections recover beta
  expect_lt(abs(run$ols_coef[["m"]] -
                  beta * sqrt(0.8 * (1 / 1.2))), 0.04)
  expect_lt(abs(run$corrected_ols_coef[["m"]] - beta), 0.05)
  expect_lt(abs(run$corrected_iv_coef[["m"]] - beta), 0.05)
})

test_that("latent outcomes (Prop 4): exact 1/sqrt(rho_Y) scaling under designs", {
  dat <- make_panel_test_data(n_units = 200L, n_periods = 8L,
                              s = 0.5, sig2U = 0.5, beta = 0.4,
                              phi = 0.7, r2_within = 0.2, seed = 515L)
  set.seed(516L)
  n <- length(dat$Y)
  yh1 <- dat$Y + rnorm(n, sd = 0.5 * stats::sd(dat$Y))
  yh2 <- dat$Y + rnorm(n, sd = 0.5 * stats::sd(dat$Y))

  run_latent <- lpmec_panel_onerun(
    Y = NULL, unit = dat$unit, time = dat$time,
    Y_split_scores = cbind(yh1, yh2),
    split_scores = list(m = cbind(dat$t1, dat$t2)),
    design = "within"
  )
  # observed-outcome run on the same constructed outcome score
  zh1 <- .lpmec_zscore(yh1)
  zh2 <- .lpmec_zscore(yh2)
  yhat <- .lpmec_zscore(rowMeans(cbind(zh1, zh2)))
  run_observed <- lpmec_panel_onerun(
    Y = yhat, unit = dat$unit, time = dat$time,
    split_scores = list(m = cbind(dat$t1, dat$t2)),
    design = "within"
  )

  expect_true(is.finite(run_latent$rho_Y))
  expect_equal(run_latent$ols_coef, run_observed$ols_coef, tolerance = 1e-12)
  outcome_factor <- 1 / sqrt(run_latent$rho_Y)
  for (field in c("corrected_ols_coef", "corrected_ols_coef_split",
                  "corrected_ols_lower", "corrected_ols_upper",
                  "corrected_iv_coef", "corrected_ols_coef_local")) {
    expect_equal(run_latent[[field]],
                 run_observed[[field]] * outcome_factor, tolerance = 1e-12)
  }
})

test_that("design-local scale (Prop 2a): bridge identity and plim recovery", {
  dat <- make_panel_test_data(n_units = 500L, n_periods = 10L,
                              s = 0.5, sig2U = 0.5, beta = 0.4,
                              phi = 0.7, r2_within = 0.2, seed = 616L)
  run <- lpmec_panel_onerun(
    Y = dat$Y, unit = dat$unit, time = dat$time,
    split_scores = list(m = cbind(dat$t1, dat$t2)),
    design = "within"
  )

  # naive local slope = pooled-scale slope times the transformed score's sd
  expect_equal(run$ols_coef_local[["m"]],
               run$ols_coef[["m"]] * run$sd_design_x[["m"]],
               tolerance = 1e-12)
  # local correction divides by sqrt of the design-scale reliability
  expect_equal(
    run$corrected_ols_coef_local_split[["m"]],
    run$ols_coef_local[["m"]] / sqrt(run$design_split_rho_score[["m"]]),
    tolerance = 1e-12
  )
  # exact local-pooled bridge (Appendix I, machine precision):
  # corrected_local = corrected_pooled * sd(L x_hat) * sqrt(rho_design/rho_pooled)
  expect_equal(
    run$corrected_ols_coef_local_split[["m"]],
    run$corrected_ols_coef_split[["m"]] * run$sd_design_x[["m"]] *
      sqrt(run$design_split_rho_score[["m"]] / run$split_rho_score[["m"]]),
    tolerance = 1e-10
  )

  # plim: the local estimand is beta * sd(design-transformed true X)
  x_within <- .lpmec_panel_transform(dat$X, unit = dat$unit,
                                     design = "within")
  expect_lt(abs(run$corrected_ols_coef_local[["m"]] -
                  dat$beta * stats::sd(x_within)), 0.05)
})
