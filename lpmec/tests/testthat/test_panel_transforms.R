# Skip all tests on CRAN to avoid timeouts
skip_on_cran()

# Unbalanced panel with time gaps and sprinkled NAs, used for brute-force
# parity checks of the design transforms.
make_unbalanced_panel <- function() {
  unit <- c(rep("A", 8), rep("B", 6), rep("C", 5))
  time <- c(1:8, c(1, 2, 4, 7, 8, 10), c(3, 5, 6, 9, 12))
  set.seed(11)
  v <- rnorm(length(unit))
  w <- rnorm(length(unit))
  v[c(2, 10, 17)] <- NA
  list(unit = unit, time = time, v = v, w = w)
}

# Independent slow implementation of NA-safe group demeaning.
brute_group_demean <- function(x, groups) {
  out <- x
  for (level in unique(groups)) {
    idx <- groups == level
    out[idx] <- x[idx] - mean(x[idx], na.rm = TRUE)
  }
  out
}

test_that("within transform matches brute force on an unbalanced panel with gaps", {
  panel <- make_unbalanced_panel()
  got <- .lpmec_panel_transform(panel$v, unit = panel$unit, time = panel$time,
                                design = "within")
  expect_equal(got, brute_group_demean(panel$v, panel$unit), tolerance = 1e-12)

  # matrix input transforms each column exactly like the vector path
  got_matrix <- .lpmec_panel_transform(cbind(v = panel$v, w = panel$w),
                                       panel$unit, panel$time,
                                       design = "within")
  expect_equal(unname(got_matrix[, "v"]), got, tolerance = 1e-12)
  expect_equal(unname(got_matrix[, "w"]),
               brute_group_demean(panel$w, panel$unit),
               tolerance = 1e-12)

  # pooled is the identity
  expect_equal(
    .lpmec_panel_transform(panel$v, panel$unit, design = "pooled"),
    panel$v
  )
})

test_that("iterated two-way demeaning matches brute force on an unbalanced panel", {
  panel <- make_unbalanced_panel()
  got <- .lpmec_panel_transform(panel$v, panel$unit, panel$time,
                                design = "twoway")
  expected <- panel$v
  for (iteration in 1:25) {
    expected <- brute_group_demean(expected, panel$unit)
    expected <- brute_group_demean(expected, panel$time)
  }
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("balanced two-way iterated demeaning equals the one-pass formula to 1e-10", {
  dat <- make_panel_test_data(n_units = 20L, n_periods = 8L)
  M <- dat$t1_mat
  one_pass <- M - rowMeans(M) - rep(colMeans(M), each = nrow(M)) + mean(M)
  iterated <- .lpmec_panel_transform(dat$t1, dat$unit, dat$time,
                                     design = "twoway")
  expect_equal(iterated, as.vector(one_pass), tolerance = 1e-10)
})

test_that("exact-gap differencing places NA where no exact t-k partner exists", {
  panel <- make_unbalanced_panel()
  b_rows <- which(panel$unit == "B")  # times 1, 2, 4, 7, 8, 10

  d1 <- .lpmec_panel_transform(panel$w, panel$unit, panel$time,
                               design = "difference", diff_k = 1L)
  expect_true(is.na(d1[b_rows[1]]))                              # t = 1
  expect_equal(d1[b_rows[2]], panel$w[b_rows[2]] - panel$w[b_rows[1]])  # t = 2
  expect_true(is.na(d1[b_rows[3]]))                              # t = 4, no t = 3
  expect_true(is.na(d1[b_rows[4]]))                              # t = 7, no t = 6
  expect_equal(d1[b_rows[5]], panel$w[b_rows[5]] - panel$w[b_rows[4]])  # t = 8
  expect_true(is.na(d1[b_rows[6]]))                              # t = 10, no t = 9

  # k = 2: t = 4 pairs with the exact t - 2 partner (t = 2) even though it is
  # not two rows back within the unit
  d2 <- .lpmec_panel_transform(panel$w, panel$unit, panel$time,
                               design = "difference", diff_k = 2L)
  expect_equal(d2[b_rows[3]], panel$w[b_rows[3]] - panel$w[b_rows[2]])
  expect_true(is.na(d2[b_rows[2]]))                              # t = 2, no t = 0

  # brute force over all rows, k = 2
  brute_d2 <- vapply(seq_along(panel$w), function(i) {
    partner <- which(panel$unit == panel$unit[i] &
                       panel$time == panel$time[i] - 2)
    if (length(partner) == 1L) panel$w[i] - panel$w[partner] else NA_real_
  }, numeric(1L))
  expect_equal(d2, brute_d2, tolerance = 1e-12)

  # brute force with NAs in the values, k = 1: NA on either side propagates
  brute_d1_v <- vapply(seq_along(panel$v), function(i) {
    partner <- which(panel$unit == panel$unit[i] &
                       panel$time == panel$time[i] - 1)
    if (length(partner) == 1L) panel$v[i] - panel$v[partner] else NA_real_
  }, numeric(1L))
  got_d1_v <- .lpmec_panel_transform(panel$v, panel$unit, panel$time,
                                     design = "difference", diff_k = 1L)
  expect_equal(got_d1_v, brute_d1_v, tolerance = 1e-12)
})

test_that("correlations with fewer than min_cor_n complete pairs return NA", {
  set.seed(3)
  a <- rnorm(25)
  b <- a + rnorm(25)
  unit <- rep(1:5, each = 5)

  gated <- .lpmec_cor_transformed(a, b, unit = unit, design = "pooled")
  expect_true(is.na(gated[["r"]]))
  expect_equal(gated[["n"]], 25)

  open <- .lpmec_cor_transformed(a, b, unit = unit, design = "pooled",
                                 min_cor_n = 10L)
  expect_equal(open[["r"]], stats::cor(a, b), tolerance = 1e-12)
  expect_equal(open[["n"]], 25)

  # the gate counts complete pairs, not rows
  a2 <- c(a, rep(NA_real_, 10))
  b2 <- c(b, rnorm(10))
  unit2 <- rep(1:7, each = 5)
  gated2 <- .lpmec_cor_transformed(a2, b2, unit = unit2, design = "pooled")
  expect_true(is.na(gated2[["r"]]))
  expect_equal(gated2[["n"]], 25)
})

test_that("cluster resampler issues fresh pseudo-ids for repeated draws", {
  unit <- rep(letters[1:6], times = c(3, 4, 2, 5, 3, 3))
  set.seed(7)
  resample <- .lpmec_resample_clusters(unit)

  expect_equal(length(resample$indices), length(resample$pseudo_unit))
  # one fresh pseudo-id per draw slot, even when a source unit repeats
  expect_equal(length(unique(resample$pseudo_unit)), 6L)
  source_unit <- sub("_[0-9]+$", "", resample$pseudo_unit)
  expect_equal(source_unit, unit[resample$indices])

  draw_sources <- sub("_[0-9]+$", "", unique(resample$pseudo_unit))
  expect_true(any(duplicated(draw_sources)))  # this seed repeats a unit

  # each pseudo-id reproduces the full row block of its source unit
  for (pseudo_id in unique(resample$pseudo_unit)) {
    rows <- resample$indices[resample$pseudo_unit == pseudo_id]
    expect_equal(rows, which(unit == sub("_[0-9]+$", "", pseudo_id)))
  }
})

test_that("panel OLS and IV slopes match their covariance forms on transformed data", {
  dat <- make_panel_test_data(n_units = 25L, n_periods = 10L)
  h1 <- as.vector(scale(dat$t1))
  h2 <- as.vector(scale(dat$t2))
  v1 <- .lpmec_panel_transform(h1, dat$unit, dat$time, design = "twoway")
  v2 <- .lpmec_panel_transform(h2, dat$unit, dat$time, design = "twoway")
  vy <- .lpmec_panel_transform(dat$Y, dat$unit, dat$time, design = "twoway")

  ols <- .lpmec_panel_ols(vy, v1, cluster = dat$unit)
  expect_equal(ols$coef, stats::cov(v1, vy) / stats::var(v1),
               tolerance = 1e-8)
  expect_true(is.finite(ols$se) && ols$se > 0)
  expect_equal(ols$n_obs, length(vy))
  expect_equal(ols$n_clusters, 25L)

  iv <- .lpmec_panel_iv(vy, v1, v2, cluster = dat$unit)
  expect_equal(iv$coef, stats::cov(v2, vy) / stats::cov(v2, v1),
               tolerance = 1e-8)
  expect_true(is.finite(iv$se) && iv$se > 0)
  expect_true(is.finite(iv$first_stage_fstat) && iv$first_stage_fstat > 0)
  expect_equal(iv$n_obs, length(vy))
})

test_that("panel input preparation enforces design requirements", {
  expect_error(
    .lpmec_prepare_panel_inputs(10L, unit = NULL, design = "within"),
    "unit"
  )
  expect_error(
    .lpmec_prepare_panel_inputs(10L, unit = rep(1:2, 5), design = "difference"),
    "time"
  )
  expect_error(
    .lpmec_prepare_panel_inputs(4L, unit = c("a", "a", "b", "b"),
                                time = c(1, 1, 1, 2), design = "within"),
    "Duplicate"
  )
  expect_error(
    .lpmec_prepare_panel_inputs(10L, unit = rep(1:2, 5), design = "banana"),
    "design"
  )

  prep <- .lpmec_prepare_panel_inputs(4L, unit = NULL, design = "pooled")
  expect_equal(prep$unit, as.character(1:4))
  expect_false(prep$unit_supplied)
})
