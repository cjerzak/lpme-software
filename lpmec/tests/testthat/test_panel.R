# Skip all tests on CRAN to avoid timeouts
skip_on_cran()

# Shared fixture: balanced panel with two parallel split scores whose
# two-way/within design reliabilities sit comfortably above the 0.05 floor.
panel_dat <- make_panel_test_data(n_units = 40L, n_periods = 10L,
                                  sig2U = 0.25, seed = 42L)

# Spearman-Brown step-up used by the split-based OLS variant
sb2 <- function(r) 2 * r / (1 + r)

test_that(".lpmec_panel_corrections implements the correction algebra exactly", {
  measure_names <- c("a", "b", "c")
  ols_coef <- c(a = 0.10, b = 0.20, c = 0.30)
  iv_coef_a <- c(a = 0.15, b = NA_real_, c = 0.25)
  iv_coef_b <- c(a = 0.17, b = NA_real_, c = 0.21)
  cross_iv_coef <- c(a_by_b = 0.11, b_by_a = 0.12, a_by_c = 0.13,
                     c_by_a = 0.14, b_by_c = 0.15, c_by_b = 0.16)
  split_pooled <- c(a = 0.8, b = NA_real_, c = 0.9)
  split_design <- c(a = 0.4, b = NA_real_, c = 0.02)  # c floored
  triad_pooled <- c(a = 0.7, b = 0.6, c = 0.5)
  triad_design <- c(a = 0.35, b = 0.30, c = 0.25)
  pair_na <- c(a = NA_real_, b = NA_real_, c = NA_real_)
  cross_pooled <- matrix(c(1, 0.6, 0.5,
                           0.6, 1, 0.4,
                           0.5, 0.4, 1), 3, 3,
                         dimnames = list(measure_names, measure_names))

  warnings_seen <- testthat::capture_warnings(
    out <- .lpmec_panel_corrections(
      ols_coef = ols_coef,
      iv_coef_a = iv_coef_a,
      iv_coef_b = iv_coef_b,
      cross_iv_coef = cross_iv_coef,
      measure_names = measure_names,
      split_pooled = split_pooled,
      split_design = split_design,
      triad_pooled = triad_pooled,
      triad_design = triad_design,
      pair_pooled = pair_na,
      pair_design = pair_na,
      cross_pooled = cross_pooled,
      min_reliability = 0.05,
      warn = TRUE
    )
  )
  expect_length(warnings_seen, 1L)
  expect_match(warnings_seen, "min_reliability")
  expect_match(warnings_seen, "c \\(split\\)")

  # corrected OLS = b * sqrt(rho_pooled) / rho_design, per variant; the
  # split variant uses the Spearman-Brown step-up of both correlations so
  # the reliabilities match the full score entered in the regression
  expect_equal(out$corrected_ols_coef_split[["a"]],
               0.10 * sqrt(sb2(0.8)) / sb2(0.4), tolerance = 1e-12)
  expect_true(is.na(out$corrected_ols_coef_split[["b"]]))  # unavailable
  expect_true(is.na(out$corrected_ols_coef_split[["c"]]))  # floored
  expect_equal(out$corrected_ols_coef_triad[["b"]], 0.20 * sqrt(0.6) / 0.30,
               tolerance = 1e-12)
  expect_equal(out$corrected_ols_coef_triad[["c"]], 0.30 * sqrt(0.5) / 0.25,
               tolerance = 1e-12)
  expect_true(all(is.na(out$corrected_ols_coef_pair)))  # M = 3: no pair

  # headline = split when available, else triad; range = [min, max]
  expect_equal(unname(out$corrected_ols_source), c("split", "triad", "triad"))
  expect_equal(out$corrected_ols_coef[["a"]],
               out$corrected_ols_coef_split[["a"]], tolerance = 1e-12)
  expect_equal(out$corrected_ols_coef[["c"]],
               out$corrected_ols_coef_triad[["c"]], tolerance = 1e-12)
  candidates_a <- c(0.10 * sqrt(sb2(0.8)) / sb2(0.4),
                    0.10 * sqrt(0.7) / 0.35)
  expect_equal(out$corrected_ols_lower[["a"]], min(candidates_a),
               tolerance = 1e-12)
  expect_equal(out$corrected_ols_upper[["a"]], max(candidates_a),
               tolerance = 1e-12)

  # IV multiplies by sqrt of the POOLED split correlation (both directions,
  # then averaged) -- never divided, never floored
  expect_equal(out$corrected_iv_coef_a[["a"]], 0.15 * sqrt(0.8),
               tolerance = 1e-12)
  expect_equal(out$corrected_iv_coef_b[["a"]], 0.17 * sqrt(0.8),
               tolerance = 1e-12)
  expect_equal(out$corrected_iv_coef_within[["a"]],
               (0.15 + 0.17) / 2 * sqrt(0.8), tolerance = 1e-12)

  # cross-measure IV at M >= 3: target-oriented (Prop 3c) -- each m_by_l
  # multiplied by sqrt of the REGRESSOR's pooled triad reliability, never
  # by sqrt of the pairwise correlation
  expect_equal(out$corrected_cross_iv_coef[["a_by_b"]], 0.11 * sqrt(0.7),
               tolerance = 1e-12)
  expect_equal(out$corrected_cross_iv_coef[["b_by_a"]], 0.12 * sqrt(0.6),
               tolerance = 1e-12)
  expect_equal(out$corrected_cross_iv_coef[["c_by_a"]], 0.14 * sqrt(0.5),
               tolerance = 1e-12)
  expect_equal(out$corrected_cross_iv_pair[["a_x_b"]],
               (0.11 * sqrt(0.7) + 0.12 * sqrt(0.6)) / 2, tolerance = 1e-12)
  expect_equal(out$corrected_iv_coef_cross[["a"]],
               mean(c(0.11 * sqrt(0.7), 0.13 * sqrt(0.7))), tolerance = 1e-12)

  # headline IV: within-split when available, else cross-measure
  expect_equal(out$corrected_iv_coef[["a"]],
               out$corrected_iv_coef_within[["a"]], tolerance = 1e-12)
  expect_equal(out$corrected_iv_coef[["b"]],
               out$corrected_iv_coef_cross[["b"]], tolerance = 1e-12)
  expect_equal(out$corrected_iv_source[["a"]], "within_split")
  expect_equal(out$corrected_iv_source[["b"]], "cross_measure")
})

test_that("two-measure parallel-pair fallback feeds the OLS correction", {
  measure_names <- c("a", "b")
  na2 <- c(a = NA_real_, b = NA_real_)
  out <- .lpmec_panel_corrections(
    ols_coef = c(a = 0.10, b = 0.20),
    iv_coef_a = na2,
    iv_coef_b = na2,
    cross_iv_coef = c(a_by_b = 0.11, b_by_a = 0.12),
    measure_names = measure_names,
    split_pooled = na2,
    split_design = na2,
    triad_pooled = na2,
    triad_design = na2,
    pair_pooled = c(a = 0.66, b = 0.66),
    pair_design = c(a = 0.33, b = 0.33),
    cross_pooled = matrix(c(1, 0.66, 0.66, 1), 2, 2,
                          dimnames = list(measure_names, measure_names)),
    min_reliability = 0.05,
    warn = TRUE
  )
  expect_equal(out$corrected_ols_coef_pair[["a"]], 0.10 * sqrt(0.66) / 0.33,
               tolerance = 1e-12)
  expect_equal(unname(out$corrected_ols_source), c("pair", "pair"))
  expect_equal(out$corrected_ols_coef, out$corrected_ols_coef_pair,
               tolerance = 1e-12)
  expect_equal(out$corrected_ols_lower, out$corrected_ols_coef_pair,
               tolerance = 1e-12)
  # cross IV corrected by sqrt of the pooled pair correlation
  expect_equal(out$corrected_cross_iv_coef[["a_by_b"]], 0.11 * sqrt(0.66),
               tolerance = 1e-12)
  expect_equal(out$corrected_iv_source[["a"]], "cross_measure")
})

test_that("scores mode returns the documented onerun structure", {
  run <- lpmec_panel_onerun(
    Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
    scores = list(s1 = panel_dat$t1, s2 = panel_dat$t2),
    design = "twoway"
  )
  expect_s3_class(run, "lpmec_panel_onerun")
  expect_equal(run$measure_names, c("s1", "s2"))
  expect_equal(run$n_measures, 2L)
  expect_equal(run$design, "twoway")
  expect_equal(run$n_obs, length(panel_dat$Y))
  expect_equal(run$n_units, panel_dat$n_units)
  expect_equal(unname(run$measure_source), c("scores", "scores"))
  expect_false(any(run$has_splits))

  expect_true(all(is.finite(run$ols_coef)))
  expect_true(all(is.finite(run$ols_se)))
  expect_equal(names(run$ols_coef), c("s1", "s2"))
  expect_equal(unname(run$ols_n_clusters), rep(panel_dat$n_units, 2L))

  # scores mode: no within-measure splits, no triads (M = 2), pair available
  expect_true(all(is.na(run$split_correlation)))
  expect_true(all(is.na(run$triad)))
  expect_true(all(is.na(run$iv_coef_a)))
  expect_true(all(is.finite(run$pair_correlation)))
  expect_true(all(is.finite(run$design_pair_correlation)))
  expect_lt(run$design_pair_correlation[["s1"]], run$pair_correlation[["s1"]])

  expect_true(all(is.finite(run$corrected_ols_coef)))
  expect_equal(unname(run$corrected_ols_source), c("pair", "pair"))
  expect_true(all(is.finite(run$cross_iv_coef)))
  expect_true(all(is.finite(run$corrected_iv_coef)))
  expect_equal(unname(run$corrected_iv_source),
               c("cross_measure", "cross_measure"))
  expect_true(all(is.finite(run$first_stage_fstat)))

  # reliability table: (pooled, twoway) x (s1, s2)
  expect_equal(nrow(run$reliability), 4L)
  expect_equal(sort(unique(run$reliability$design)),
               sort(c("pooled", "twoway")))
  expect_equal(names(run$cor_matrices), c("pooled", "twoway"))

  # correction algebra ties out against the returned reliabilities
  expect_equal(
    run$corrected_ols_coef[["s1"]],
    run$ols_coef[["s1"]] * sqrt(run$pair_correlation[["s1"]]) /
      run$design_pair_correlation[["s1"]],
    tolerance = 1e-12
  )
  expect_equal(
    run$corrected_iv_coef[["s1"]],
    run$cross_iv_coef[["s1_by_s2"]] * sqrt(run$pair_correlation[["s1"]]),
    tolerance = 1e-12
  )
})

test_that("hybrid observables + scores + split_scores inputs merge (M = 3)", {
  set.seed(88)
  n <- length(panel_dat$Y)
  items <- panel_dat$X + matrix(rnorm(n * 6L, sd = 0.5), n, 6L)
  colnames(items) <- paste0("item", 1:6)
  s_score <- panel_dat$X + rnorm(n, sd = 0.5)

  run <- lpmec_panel_onerun(
    Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
    observables = list(k = items),
    scores = list(s = s_score),
    split_scores = list(m = cbind(panel_dat$t1, panel_dat$t2)),
    design = "within",
    estimation_method = "averaging"
  )
  expect_equal(run$measure_names, c("k", "s", "m"))
  expect_equal(unname(run$measure_source),
               c("observables", "scores", "split_scores"))
  expect_equal(unname(run$has_splits), c(TRUE, FALSE, TRUE))

  # triads available for all measures; splits only where halves exist
  expect_true(all(is.finite(run$triad)))
  expect_true(all(is.finite(run$design_triad)))
  expect_true(is.finite(run$split_correlation[["k"]]))
  expect_true(is.finite(run$split_correlation[["m"]]))
  expect_true(is.na(run$split_correlation[["s"]]))
  expect_true(all(is.na(run$pair_correlation)))  # M = 3: pair not used

  expect_equal(unname(run$corrected_ols_source), c("split", "triad", "split"))
  expect_equal(unname(run$corrected_iv_source),
               c("within_split", "cross_measure", "within_split"))
  expect_true(all(is.finite(run$corrected_ols_coef)))
  expect_true(all(run$corrected_ols_lower <= run$corrected_ols_upper))
  expect_true(all(is.finite(run$iv_coef[c("k", "m")])))
  expect_length(run$cross_iv_coef, 6L)  # 3 x 2 ordered pairs
  expect_length(run$corrected_cross_iv_pair, 3L)
  expect_true(all(is.finite(run$corrected_iv_coef)))
})

test_that("items mode and split_scores mode agree on all split-based output", {
  set.seed(77)
  n <- length(panel_dat$Y)
  items <- panel_dat$X + matrix(rnorm(n * 6L, sd = 0.6), n, 6L)
  colnames(items) <- paste0("item", 1:6)
  partition <- list(split1_names = c("item1", "item2", "item3"),
                    split2_names = c("item4", "item5", "item6"))

  scalar_run <- lpmec_onerun(
    Y = panel_dat$Y, observables = items,
    estimation_method = "averaging", partition = partition
  )
  run_items <- lpmec_panel_onerun(
    Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
    observables = list(k = items),
    partition = list(k = partition),
    design = "within",
    estimation_method = "averaging"
  )
  run_splits <- lpmec_panel_onerun(
    Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
    split_scores = list(k = cbind(as.numeric(scalar_run$x_est1),
                                  as.numeric(scalar_run$x_est2))),
    design = "within"
  )

  # identical half scores after pooled z-scoring and alignment
  expect_equal(unname(run_items$x_est1[, "k"]),
               unname(run_splits$x_est1[, "k"]), tolerance = 1e-10)
  expect_equal(unname(run_items$x_est2[, "k"]),
               unname(run_splits$x_est2[, "k"]), tolerance = 1e-10)

  # all split-based quantities agree across the two input modes
  for (field in c("split_correlation", "design_split_correlation",
                  "iv_coef_a", "iv_coef_b", "iv_coef",
                  "corrected_iv_coef_a", "corrected_iv_coef_b",
                  "corrected_iv_coef_within", "corrected_iv_coef")) {
    expect_equal(run_items[[field]], run_splits[[field]], tolerance = 1e-10)
  }
})

test_that("covariates are design-transformed and enter the OLS stage", {
  set.seed(99)
  n <- length(panel_dat$Y)
  c1 <- 0.5 * panel_dat$X + rnorm(n)
  run_cov <- lpmec_panel_onerun(
    Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
    scores = list(s1 = panel_dat$t1, s2 = panel_dat$t2),
    covariates = data.frame(c1 = c1),
    design = "within"
  )
  expect_equal(run_cov$covariate_names, "c1")

  # manual within-transformed regression with the covariate
  z1 <- .lpmec_zscore(panel_dat$t1)
  x_w <- .lpmec_panel_transform(z1, unit = panel_dat$unit, design = "within")
  y_w <- .lpmec_panel_transform(panel_dat$Y, unit = panel_dat$unit,
                                design = "within")
  c_w <- .lpmec_panel_transform(c1, unit = panel_dat$unit, design = "within")
  manual_fit <- stats::lm(y_w ~ x_w + c_w)
  expect_equal(unname(run_cov$ols_coef[["s1"]]),
               unname(stats::coef(manual_fit)[["x_w"]]), tolerance = 1e-10)

  # covariate changes the naive fit relative to the no-covariate run
  run_nocov <- lpmec_panel_onerun(
    Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
    scores = list(s1 = panel_dat$t1, s2 = panel_dat$t2),
    design = "within"
  )
  expect_false(isTRUE(all.equal(run_cov$ols_coef[["s1"]],
                                run_nocov$ols_coef[["s1"]])))
  expect_true(all(is.finite(run_cov$corrected_ols_coef)))
  expect_true(all(is.finite(run_cov$corrected_iv_coef)))
})

test_that("validation errors identify malformed panel inputs", {
  n <- length(panel_dat$Y)
  two_scores <- list(s1 = panel_dat$t1, s2 = panel_dat$t2)

  expect_error(
    lpmec_panel_onerun(Y = panel_dat$Y, scores = two_scores),
    "'unit' is required"
  )
  expect_error(
    lpmec_panel_onerun(Y = panel_dat$Y, unit = panel_dat$unit,
                       time = rep(1, n), scores = two_scores,
                       design = "twoway"),
    "Duplicate \\(unit, time\\)"
  )
  expect_error(
    lpmec_panel_onerun(Y = panel_dat$Y, unit = panel_dat$unit,
                       scores = two_scores, design = "difference"),
    "'time' is required"
  )
  expect_error(
    lpmec_panel_onerun(Y = panel_dat$Y, unit = panel_dat$unit,
                       scores = list(only = panel_dat$t1)),
    "single measure without split halves"
  )
  expect_error(
    lpmec_panel_onerun(Y = panel_dat$Y[-1], unit = panel_dat$unit[-1],
                       scores = two_scores),
    "one value per row"
  )
  expect_error(
    lpmec_panel_onerun(Y = panel_dat$Y, unit = panel_dat$unit,
                       scores = two_scores, design = "banana"),
    "'arg'"
  )
  expect_error(
    lpmec_panel(Y = panel_dat$Y, unit = panel_dat$unit,
                scores = two_scores, n_boot = -1L),
    "'n_boot'"
  )
  expect_error(
    lpmec_panel(Y = panel_dat$Y, unit = panel_dat$unit,
                scores = two_scores, n_partition = 0L),
    "'n_partition'"
  )
  expect_error(
    lpmec_panel(Y = panel_dat$Y, unit = panel_dat$unit[-1],
                scores = two_scores),
    "'unit' must have the same length"
  )
})

test_that("lpmec_panel coerces n_partition to 1 without items mode", {
  expect_message(
    suppressWarnings(lpmec_panel(
      Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
      scores = list(s1 = panel_dat$t1, s2 = panel_dat$t2),
      design = "twoway", n_boot = 0L, n_partition = 10L
    )),
    "n_partition"
  )
})

test_that("cluster bootstrap output has the documented shapes (n_boot = 2)", {
  scores_input <- list(s1 = panel_dat$t1, s2 = panel_dat$t2)
  agg <- suppressMessages(lpmec_panel(
    Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
    scores = scores_input, design = "twoway",
    n_boot = 2L, seed = 99
  ))
  expect_s3_class(agg, "lpmec_panel")
  expect_equal(agg$n_boot, 2L)
  expect_equal(agg$n_partition, 1L)  # coerced: no items mode
  expect_equal(agg$boot_n_failed, 0L)
  expect_equal(agg$measure_names, c("s1", "s2"))

  # point estimates equal the single-partition original run
  run <- lpmec_panel_onerun(
    Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
    scores = scores_input, design = "twoway"
  )
  expect_equal(agg$ols_coef, run$ols_coef, tolerance = 1e-12)
  expect_equal(agg$corrected_ols_coef, run$corrected_ols_coef,
               tolerance = 1e-12)
  expect_equal(agg$ols_cluster_se, run$ols_se, tolerance = 1e-12)

  # bootstrap summaries exist for every aggregated field
  for (field in c("ols_coef", "corrected_ols_coef", "corrected_iv_coef",
                  "pair_correlation", "design_pair_correlation",
                  "cross_iv_coef", "corrected_cross_iv_pair")) {
    for (suffix in c("_se", "_lower", "_upper")) {
      expect_true(paste0(field, suffix) %in% names(agg))
    }
  }
  expect_true(all(is.finite(agg$ols_coef_se)))
  expect_true(all(agg$ols_coef_lower <= agg$ols_coef_upper))

  # intermediaries: (n_boot + 1) x n_partition runs, row 1 = original
  expect_equal(agg$Intermediary_BootIndex, c(1L, 2L, 3L))
  expect_equal(agg$Intermediary_PartitionIndex, rep(1L, 3L))
  expect_equal(dim(agg$Intermediary_ols_coef), c(3L, 2L))
  expect_equal(unname(agg$Intermediary_ols_coef[1L, ]), unname(run$ols_coef),
               tolerance = 1e-12)

  # reproducible under the same seed
  agg2 <- suppressMessages(lpmec_panel(
    Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
    scores = scores_input, design = "twoway",
    n_boot = 2L, seed = 99
  ))
  expect_equal(agg2$ols_coef_se, agg$ols_coef_se, tolerance = 1e-12)
  expect_equal(agg2$Intermediary_ols_coef, agg$Intermediary_ols_coef,
               tolerance = 1e-12)

  # n_boot = 0: no bootstrap summaries
  agg0 <- suppressMessages(lpmec_panel(
    Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
    scores = scores_input, design = "twoway", n_boot = 0L
  ))
  expect_true(all(is.na(agg0$ols_coef_se)))
  expect_equal(agg0$ols_coef, run$ols_coef, tolerance = 1e-12)

  # return_intermediaries = FALSE drops the per-run matrices
  agg_lean <- suppressMessages(lpmec_panel(
    Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
    scores = scores_input, design = "twoway", n_boot = 0L,
    return_intermediaries = FALSE
  ))
  expect_false(any(grepl("^Intermediary_", names(agg_lean))))
})

test_that("items mode keeps n_partition and varies partitions across runs", {
  dat <- make_panel_test_data(n_units = 40L, n_periods = 10L, seed = 21L)
  set.seed(61)
  n <- length(dat$Y)
  items <- dat$X + matrix(rnorm(n * 6L, sd = 0.9), n, 6L)
  colnames(items) <- paste0("it", 1:6)

  agg <- suppressMessages(lpmec_panel(
    Y = dat$Y, unit = dat$unit, time = dat$time,
    observables = list(k = items),
    scores = list(s = dat$t1),
    design = "within",
    estimation_method = "averaging",
    n_boot = 1L, n_partition = 2L, seed = 3
  ))
  expect_equal(agg$n_partition, 2L)
  expect_equal(agg$Intermediary_BootIndex, c(1L, 1L, 2L, 2L))
  expect_equal(agg$Intermediary_PartitionIndex, c(1L, 2L, 1L, 2L))
  expect_equal(dim(agg$Intermediary_ols_coef), c(4L, 2L))
  # fresh random partitions change the split correlation of the item measure
  split_k <- agg$Intermediary_split_correlation[, "k"]
  expect_true(all(is.finite(split_k)))
  expect_gt(stats::sd(split_k[1:2]), 0)
})

test_that("single split-scores measure aggregates through the 1-column path", {
  dat <- make_panel_test_data(n_units = 40L, n_periods = 10L, seed = 21L)
  agg <- suppressMessages(lpmec_panel(
    Y = dat$Y, unit = dat$unit, time = dat$time,
    split_scores = list(m = cbind(dat$t1, dat$t2)),
    design = "within", n_boot = 2L, seed = 5
  ))
  expect_equal(agg$measure_names, "m")
  expect_equal(names(agg$ols_coef), "m")
  expect_true(is.finite(agg$ols_coef[["m"]]))
  expect_true(is.finite(agg$ols_coef_se[["m"]]))
  expect_true(is.finite(agg$corrected_ols_coef[["m"]]))
  expect_equal(dim(agg$Intermediary_ols_coef), c(3L, 1L))
  # no cross-measure quantities with a single measure
  expect_length(agg$cross_iv_coef, 0L)
  expect_length(agg$corrected_cross_iv_pair, 0L)
})
