# Skip all tests on CRAN to avoid timeouts
skip_on_cran()

# ---------------------------------------------------------------------------
# Acceptance gate: parity with the V2 reference estimator pipeline.
# estimate_core() below is a line-for-line re-implementation of the verified
# reference in LatentMeasures/V2/Code/sim_paper/00_common.R:
#   naive : two-way FE OLS of Y on the pooled-standardized split-1 score
#   corOLS: Prop 3a, naive * sqrt(rho_pooled_hat) / rho_within_hat
#   corIV : Prop 3b, within split-IV (instrument = split 2) * sqrt(rho_pooled)
# ---------------------------------------------------------------------------
ref_std <- function(x) (x - mean(x)) / sd(x)
ref_std_mat <- function(M) matrix(ref_std(as.vector(M)), nrow(M))
ref_tw_mat <- function(M) {
  M - rowMeans(M) - rep(colMeans(M), each = nrow(M)) + mean(M)
}
ref_estimate_core <- function(Y, t1, t2) {
  h1 <- ref_std_mat(t1)
  h2 <- ref_std_mat(t2)
  v1 <- as.vector(ref_tw_mat(h1))
  v2 <- as.vector(ref_tw_mat(h2))
  vy <- as.vector(ref_tw_mat(Y))
  rho_p <- cor(as.vector(h1), as.vector(h2))
  rho_w <- cor(v1, v2)
  naive <- cov(v1, vy) / var(v1)
  iv <- cov(v2, vy) / cov(v2, v1)
  c(naive = naive,
    corOLS = naive * sqrt(rho_p) / rho_w,
    corIV = iv * sqrt(rho_p),
    pooled = cov(as.vector(h1), as.vector(Y)) / var(as.vector(h1)),
    rho_p = rho_p,
    rho_w = rho_w)
}

test_that("parity gate: onerun matches V2 estimate_core to 1e-8 (balanced twoway)", {
  dat <- make_panel_test_data(n_units = 40L, n_periods = 10L,
                              sig2U = 0.25, seed = 42L)
  reference <- ref_estimate_core(dat$Y_mat, dat$t1_mat, dat$t2_mat)

  run <- lpmec_panel_onerun(
    Y = dat$Y, unit = dat$unit, time = dat$time,
    scores = list(s1 = dat$t1, s2 = dat$t2),
    design = "twoway"
  )

  expect_equal(run$ols_coef[["s1"]], reference[["naive"]], tolerance = 1e-8)
  expect_equal(run$corrected_ols_coef[["s1"]], reference[["corOLS"]],
               tolerance = 1e-8)
  expect_equal(run$corrected_cross_iv_coef[["s1_by_s2"]],
               reference[["corIV"]], tolerance = 1e-8)
  expect_equal(run$corrected_iv_coef[["s1"]], reference[["corIV"]],
               tolerance = 1e-8)
  expect_equal(run$pair_correlation[["s1"]], reference[["rho_p"]],
               tolerance = 1e-8)
  expect_equal(run$design_pair_correlation[["s1"]], reference[["rho_w"]],
               tolerance = 1e-8)

  # explicit max-abs-deviation report for the acceptance record
  deviations <- c(
    naive = run$ols_coef[["s1"]] - reference[["naive"]],
    corOLS = run$corrected_ols_coef[["s1"]] - reference[["corOLS"]],
    corIV = run$corrected_iv_coef[["s1"]] - reference[["corIV"]],
    rho_p = run$pair_correlation[["s1"]] - reference[["rho_p"]],
    rho_w = run$design_pair_correlation[["s1"]] - reference[["rho_w"]]
  )
  expect_lt(max(abs(deviations)), 1e-8)
})

test_that("first-difference design matches the covariance-form algebra to 1e-8", {
  dat <- make_panel_test_data(n_units = 40L, n_periods = 10L,
                              sig2U = 0.25, seed = 42L)
  h1 <- ref_std(dat$t1)
  h2 <- ref_std(dat$t2)
  row_key <- paste(dat$unit, dat$time)
  partner <- match(paste(dat$unit, dat$time - 1L), row_key)
  first_diff <- function(v) v - v[partner]
  dh1 <- first_diff(h1)
  dh2 <- first_diff(h2)
  dy <- first_diff(dat$Y)
  complete <- stats::complete.cases(dh1, dh2, dy)

  naive_ref <- cov(dh1[complete], dy[complete]) / var(dh1[complete])
  rho_p_ref <- cor(h1, h2)
  rho_d_ref <- cor(dh1[complete], dh2[complete])
  corOLS_ref <- naive_ref * sqrt(rho_p_ref) / rho_d_ref
  corIV_ref <- cov(dh2[complete], dy[complete]) /
    cov(dh2[complete], dh1[complete]) * sqrt(rho_p_ref)

  run <- lpmec_panel_onerun(
    Y = dat$Y, unit = dat$unit, time = dat$time,
    scores = list(s1 = dat$t1, s2 = dat$t2),
    design = "difference", diff_k = 1L
  )

  expect_equal(run$ols_coef[["s1"]], naive_ref, tolerance = 1e-8)
  expect_equal(run$corrected_ols_coef[["s1"]], corOLS_ref, tolerance = 1e-8)
  expect_equal(run$corrected_cross_iv_coef[["s1_by_s2"]], corIV_ref,
               tolerance = 1e-8)
  expect_equal(run$pair_correlation[["s1"]], rho_p_ref, tolerance = 1e-8)
  expect_equal(run$design_pair_correlation[["s1"]], rho_d_ref,
               tolerance = 1e-8)
  # exactly the first period lacks a t - 1 partner on a balanced panel
  expect_equal(unname(run$ols_n_obs[["s1"]]),
               dat$n_units * (dat$n_periods - 1L))
})

test_that("plim recovery at N = 600, T = 8, sig2U = 0.1 (finite-T theory)", {
  beta <- 0.4
  phi <- 0.95
  s <- 0.5
  sig2U <- 0.1
  n_periods <- 8L

  # finite-T theory re-derived here (V2 sim_paper/00_common.R conventions)
  ar1_within_var <- function(phi, v, T_periods) {
    g <- v * phi^abs(outer(seq_len(T_periods), seq_len(T_periods), "-"))
    v - mean(g)
  }
  rho_fe_exact <- function(s, sig2U, T_periods) {
    Vw <- ar1_within_var(phi, s, T_periods)
    Vw / (Vw + sig2U * (T_periods - 1L) / T_periods)
  }
  rho_fe <- rho_fe_exact(s, sig2U, n_periods)
  rho_pooled <- 1 / (1 + sig2U)
  naive_plim <- beta * sqrt(1 + sig2U) * rho_fe

  dat <- make_panel_test_data(n_units = 600L, n_periods = n_periods,
                              s = s, sig2U = sig2U, beta = beta,
                              phi = phi, seed = 11L)
  run <- lpmec_panel_onerun(
    Y = dat$Y, unit = dat$unit, time = dat$time,
    scores = list(s1 = dat$t1, s2 = dat$t2),
    design = "twoway"
  )

  # observed deviations (seed 11): rho_p 0.0002, rho_w 0.0005, naive 0.008,
  # corOLS 0.017, corIV 0.005 -- tolerances leave a ~3x margin
  expect_equal(run$pair_correlation[["s1"]], rho_pooled, tolerance = 0.02)
  expect_equal(run$design_pair_correlation[["s1"]], rho_fe, tolerance = 0.05)
  expect_lt(abs(run$ols_coef[["s1"]] - naive_plim), 0.03)
  expect_lt(abs(run$corrected_ols_coef[["s1"]] - beta), 0.05)
  expect_lt(abs(run$corrected_iv_coef[["s1"]] - beta), 0.05)

  # the FE design amplifies attenuation: naive under FE is far below the
  # pooled-design plim beta / sqrt(1 + sig2U)
  expect_lt(run$ols_coef[["s1"]], beta / sqrt(1 + sig2U) - 0.1)
})

test_that("outputs are invariant to flipping the sign of one input measure", {
  dat <- make_panel_test_data(n_units = 40L, n_periods = 10L,
                              sig2U = 0.25, seed = 42L)
  base <- lpmec_panel_onerun(
    Y = dat$Y, unit = dat$unit, time = dat$time,
    scores = list(s1 = dat$t1, s2 = dat$t2),
    design = "twoway"
  )
  flipped <- lpmec_panel_onerun(
    Y = dat$Y, unit = dat$unit, time = dat$time,
    scores = list(s1 = dat$t1, s2 = -dat$t2),
    design = "twoway"
  )
  expect_true(flipped$sign_flipped[["s2"]])
  expect_false(base$sign_flipped[["s2"]])
  for (field in c("ols_coef", "corrected_ols_coef", "pair_correlation",
                  "design_pair_correlation", "cross_iv_coef",
                  "corrected_cross_iv_coef", "corrected_iv_coef")) {
    expect_equal(flipped[[field]], base[[field]], tolerance = 1e-10)
  }
})

test_that("design reliability below the floor yields NA plus one warning", {
  set.seed(33)
  n_units <- 50L
  n_periods <- 10L
  unit <- rep(seq_len(n_units), times = n_periods)
  time <- rep(seq_len(n_periods), each = n_units)
  # halves share only a unit-level component: pooled reliability is healthy,
  # within-design reliability is ~0 (below the 0.05 floor)
  unit_component <- rnorm(n_units)
  half1 <- unit_component[unit] + rnorm(length(unit))
  half2 <- unit_component[unit] + rnorm(length(unit))
  Y <- rnorm(length(unit))

  warnings_seen <- testthat::capture_warnings(
    run <- lpmec_panel_onerun(
      Y = Y, unit = unit, time = time,
      split_scores = list(m = cbind(half1, half2)),
      design = "within"
    )
  )
  expect_length(warnings_seen, 1L)
  expect_match(warnings_seen, "min_reliability")

  expect_gt(run$split_correlation[["m"]], 0.05)     # pooled: identified
  expect_lt(run$design_split_correlation[["m"]], 0.05)  # design: floored
  expect_true(is.na(run$corrected_ols_coef[["m"]]))
  expect_true(is.na(run$corrected_ols_coef_split[["m"]]))
  expect_true(is.na(run$corrected_ols_lower[["m"]]))
  expect_true(is.na(run$corrected_ols_upper[["m"]]))
  # IV multiplies by sqrt(rho_pooled) and stays defined where OLS divides
  # by a floored design reliability
  expect_true(is.finite(run$corrected_iv_coef[["m"]]))
})
