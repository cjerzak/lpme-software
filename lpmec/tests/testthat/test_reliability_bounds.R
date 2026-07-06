# Skip all tests on CRAN to avoid timeouts
skip_on_cran()

# Prop-6 mini-sim worlds. Measure 1 is built from two half scores that share
# a within-measure error u (so its split correlation overstates the score
# reliability); measures 2 and 3 are full scores. In the divergence world,
# measures 2 and 3 also share a cross-measure method error, which inflates
# r_23 and biases measure 1's triad downward.
#
# Population values (Var(X) = 1):
#   measure-1 score ~ X + u + (w1 + w2) / 2, so rho*_1 = 1 / (1 + 0.3 + 0.2)
#   split_1 = (1 + 0.3) / (1 + 0.3 + 0.4)
simulate_prop6_world <- function(n, shared_cross_error = FALSE, seed = 21) {
  set.seed(seed)
  X <- rnorm(n)
  u <- rnorm(n, sd = sqrt(0.3))
  h1 <- X + u + rnorm(n, sd = sqrt(0.4))
  h2 <- X + u + rnorm(n, sd = sqrt(0.4))
  if (shared_cross_error) {
    c_shared <- rnorm(n, sd = sqrt(0.25))
    s2 <- X + c_shared + rnorm(n, sd = sqrt(0.25))
    s3 <- X + c_shared + rnorm(n, sd = sqrt(0.25))
    rho_star_2 <- 1 / 1.5
  } else {
    s2 <- X + rnorm(n, sd = sqrt(0.25))
    s3 <- X + rnorm(n, sd = sqrt(0.5))
    rho_star_2 <- 1 / 1.25
  }
  list(
    X = X, h1 = h1, h2 = h2, s2 = s2, s3 = s3,
    rho_star_1 = 1 / 1.5,
    rho_star_2 = rho_star_2,
    split_plim_1 = 1.3 / 1.7
  )
}

test_that("spearman-brown step-up and triad identities are exact", {
  expect_equal(.lpmec_spearman_brown(0.342, 2), 2 * 0.342 / (1 + 0.342),
               tolerance = 1e-12)
  expect_equal(.lpmec_spearman_brown(0.5, 4), 4 * 0.5 / (1 + 3 * 0.5),
               tolerance = 1e-12)
  expect_equal(.lpmec_spearman_brown(c(0.2, 0.8), 2),
               c(2 * 0.2 / 1.2, 2 * 0.8 / 1.8), tolerance = 1e-12)
  expect_true(is.na(.lpmec_spearman_brown(NA_real_, 2)))
  expect_error(.lpmec_spearman_brown(0.5, 0), "factor")

  R <- matrix(c(1, 0.6, 0.5,
                0.6, 1, 0.4,
                0.5, 0.4, 1), 3, 3,
              dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  triads <- .lpmec_triad_reliabilities(R)
  expect_equal(unname(triads["a"]), 0.6 * 0.5 / 0.4, tolerance = 1e-12)
  expect_equal(unname(triads["b"]), 0.6 * 0.4 / 0.5, tolerance = 1e-12)
  expect_equal(unname(triads["c"]), 0.5 * 0.4 / 0.6, tolerance = 1e-12)

  # fewer than 3 measures: all NA
  expect_true(all(is.na(.lpmec_triad_reliabilities(R[1:2, 1:2]))))
})

test_that("Prop-6 mini-sim: triads recover rho* and splits sit above it", {
  world <- simulate_prop6_world(4000, shared_cross_error = FALSE, seed = 21)
  bounds <- lpmec_reliability_bounds(
    split_scores = list(m1 = cbind(world$h1, world$h2)),
    scores = list(m2 = world$s2, m3 = world$s3)
  )
  rel <- bounds$reliability
  row1 <- rel[rel$measure == "m1", ]
  row2 <- rel[rel$measure == "m2", ]

  # triad consistency under independent cross-measure errors (tol 0.03)
  expect_equal(row1$triad, world$rho_star_1, tolerance = 0.03)
  expect_equal(row2$triad, world$rho_star_2, tolerance = 0.03)

  # shared within-measure half error: split is an upper bound of rho*
  expect_equal(row1$split_correlation, world$split_plim_1, tolerance = 0.03)
  expect_gt(row1$split_correlation, world$rho_star_1)

  # bounds are the [min, max] of the finite candidates
  expect_equal(row1$rho_lo, min(row1$triad, row1$split_correlation))
  expect_equal(row1$rho_hi, max(row1$triad, row1$split_correlation))
})

test_that("Prop-6 mini-sim: divergence world orders triad < rho* < split", {
  world <- simulate_prop6_world(4000, shared_cross_error = TRUE, seed = 22)
  bounds <- lpmec_reliability_bounds(
    split_scores = list(m1 = cbind(world$h1, world$h2)),
    scores = list(m2 = world$s2, m3 = world$s3)
  )
  row1 <- bounds$reliability[bounds$reliability$measure == "m1", ]

  # population values: triad_1 = 8/15, rho*_1 = 2/3, split_1 = 13/17
  expect_lt(row1$triad, world$rho_star_1 - 0.05)
  expect_gt(row1$split_correlation, world$rho_star_1 + 0.05)
  expect_lt(row1$rho_lo, row1$rho_hi)
  expect_equal(row1$rho_lo, row1$triad)
  expect_equal(row1$rho_hi, row1$split_correlation)
})

test_that("two measures yield NA triads and split-only bounds", {
  set.seed(31)
  n <- 300
  X <- rnorm(n)
  splits <- list(
    a = cbind(X + rnorm(n, sd = 0.5), X + rnorm(n, sd = 0.5)),
    b = cbind(X + rnorm(n, sd = 0.7), X + rnorm(n, sd = 0.7))
  )
  bounds <- lpmec_reliability_bounds(split_scores = splits)
  rel <- bounds$reliability
  expect_true(all(is.na(rel$triad)))
  expect_true(all(is.finite(rel$split_correlation)))
  expect_equal(rel$rho_lo, rel$split_correlation)
  expect_equal(rel$rho_hi, rel$split_correlation)
})

test_that("multiple designs produce one reliability block per design", {
  dat <- make_panel_test_data(n_units = 30L, n_periods = 20L)
  set.seed(41)
  s2 <- dat$X + rnorm(length(dat$X), sd = sqrt(dat$sig2U))
  s3 <- dat$X + rnorm(length(dat$X), sd = sqrt(dat$sig2U))
  designs <- c("pooled", "within", "twoway", "difference")

  # the difference design attenuates reliabilities toward (possibly below)
  # the floor, so the single floor warning may legitimately fire here
  bounds <- suppressWarnings(lpmec_reliability_bounds(
    split_scores = list(m1 = cbind(dat$t1, dat$t2)),
    scores = list(m2 = s2, m3 = s3),
    unit = dat$unit,
    time = dat$time,
    designs = designs,
    diff_k = 1L
  ))
  rel <- bounds$reliability

  expect_equal(nrow(rel), 4L * 3L)
  expect_equal(unique(rel$design), designs)
  expect_equal(names(bounds$cor_matrices), designs)
  expect_equal(dim(bounds$cor_matrices[["within"]]), c(3L, 3L))

  m1 <- rel[rel$measure == "m1", ]
  expect_true(all(is.finite(m1$split_correlation)))
  expect_true(all(is.finite(m1$triad)))
  # design reliabilities attenuate relative to pooled under FE transforms
  expect_gt(m1$split_correlation[m1$design == "pooled"],
            m1$split_correlation[m1$design == "twoway"])
  expect_gt(m1$split_correlation[m1$design == "pooled"],
            m1$split_correlation[m1$design == "difference"])
})

test_that("bootstrap adds uncertainty columns and intermediaries", {
  set.seed(51)
  n <- 200
  X <- rnorm(n)
  scores <- list(
    a = X + rnorm(n, sd = 0.5),
    b = X + rnorm(n, sd = 0.5),
    c = X + rnorm(n, sd = 0.5)
  )
  bounds <- lpmec_reliability_bounds(scores = scores, n_boot = 4L, seed = 99)
  rel <- bounds$reliability

  for (field in c("split_correlation", "triad", "rho_lo", "rho_hi")) {
    for (suffix in c("_se", "_lower", "_upper")) {
      expect_true(paste0(field, suffix) %in% names(rel))
    }
  }
  for (field in c("triad", "rho_lo", "rho_hi")) {
    expect_true(all(is.finite(rel[[paste0(field, "_se")]])))
    expect_true(all(rel[[paste0(field, "_lower")]] <=
                      rel[[paste0(field, "_upper")]]))
  }

  expect_equal(bounds$n_boot, 4L)
  expect_equal(bounds$n_boot_failed, 0L)
  expect_equal(bounds$Intermediary_BootIndex, seq_len(5L))
  expect_equal(dim(bounds$Intermediary_triad), c(5L, 3L))
  expect_equal(unname(bounds$Intermediary_triad[1L, ]), rel$triad)

  # reproducible under the same seed
  bounds2 <- lpmec_reliability_bounds(scores = scores, n_boot = 4L, seed = 99)
  expect_equal(bounds2$reliability, rel)
})

test_that("observables-mode measures merge with scores and expose split scores", {
  set.seed(61)
  n <- 150
  X <- rnorm(n)
  items <- X + matrix(rnorm(n * 6, sd = 0.8), n, 6)
  colnames(items) <- paste0("item", 1:6)
  s2 <- X + rnorm(n, sd = 0.6)

  bounds <- lpmec_reliability_bounds(
    observables = list(k = items),
    scores = list(s = s2),
    estimation_method = "averaging",
    seed = 5
  )
  expect_s3_class(bounds, "lpmec_reliability_bounds")
  expect_equal(bounds$measure_names, c("k", "s"))
  expect_equal(unname(bounds$measure_source),
               c("observables", "scores"))
  expect_true(all(is.finite(bounds$x_est1[, "k"])))
  expect_true(all(is.na(bounds$x_est1[, "s"])))

  rel <- bounds$reliability
  expect_true(is.finite(rel$split_correlation[rel$measure == "k"]))
  expect_true(is.na(rel$split_correlation[rel$measure == "s"]))
  expect_true(all(is.na(rel$triad)))  # only 2 measures
})

test_that("sign alignment makes cross-measure reliabilities orientation-invariant", {
  world <- simulate_prop6_world(2000, shared_cross_error = FALSE, seed = 71)
  base <- lpmec_reliability_bounds(
    split_scores = list(m1 = cbind(world$h1, world$h2)),
    scores = list(m2 = world$s2, m3 = world$s3)
  )
  # negate one input measure: reliabilities and bounds must be unchanged
  flipped <- lpmec_reliability_bounds(
    split_scores = list(m1 = cbind(world$h1, world$h2)),
    scores = list(m2 = -world$s2, m3 = world$s3)
  )
  expect_equal(flipped$reliability$triad, base$reliability$triad,
               tolerance = 1e-8)
  expect_equal(flipped$reliability$split_correlation,
               base$reliability$split_correlation, tolerance = 1e-8)
  expect_equal(flipped$reliability$rho_lo, base$reliability$rho_lo,
               tolerance = 1e-8)
  expect_equal(flipped$reliability$rho_hi, base$reliability$rho_hi,
               tolerance = 1e-8)
  expect_equal(flipped$cor_matrices$pooled, base$cor_matrices$pooled,
               tolerance = 1e-8)

  # no flips are needed for consistently oriented inputs; the negated run
  # flips at least one column and ends internally consistent (all pairwise
  # correlations positive), as does the base run
  expect_false(any(base$sign_flipped))
  expect_true(any(flipped$sign_flipped))
  for (result in list(base, flipped)) {
    aligned_cor <- stats::cor(result$x_est, use = "pairwise.complete.obs")
    expect_true(all(aligned_cor > 0))
  }

  # halves are aligned with their own measure score
  expect_gt(stats::cor(flipped$x_est1[, "m1"], flipped$x_est[, "m1"],
                       use = "pairwise.complete.obs"), 0)
})

test_that("reliabilities below the floor are dropped from bounds with one warning", {
  tab <- data.frame(
    design = c("within", "within"),
    measure = c("a", "b"),
    split_correlation = c(0.02, 0.6),
    triad = c(0.4, -0.1),
    stringsAsFactors = FALSE
  )
  warnings_seen <- testthat::capture_warnings(
    out <- .lpmec_bounds_from_reliability(tab, min_reliability = 0.05)
  )
  expect_length(warnings_seen, 1L)
  expect_match(warnings_seen, "min_reliability")
  expect_equal(out$rho_lo[1], 0.4)
  expect_equal(out$rho_hi[1], 0.4)
  expect_equal(out$rho_lo[2], 0.6)
  expect_equal(out$rho_hi[2], 0.6)

  # no finite candidate above the floor: bounds are NA
  tab_na <- data.frame(
    design = "within", measure = "c",
    split_correlation = 0.01, triad = NA_real_,
    stringsAsFactors = FALSE
  )
  suppressWarnings(out_na <- .lpmec_bounds_from_reliability(tab_na))
  expect_true(is.na(out_na$rho_lo))
  expect_true(is.na(out_na$rho_hi))
})

test_that("input validation catches duplicate names and missing requirements", {
  expect_error(lpmec_reliability_bounds(), "At least one")

  set.seed(81)
  n <- 120
  X <- rnorm(n)
  s <- X + rnorm(n)
  # Same name in 'scores' and 'split_scores' merges (full score + halves).
  s_h1 <- X + rnorm(n)
  s_h2 <- X + rnorm(n)
  merged <- lpmec_reliability_bounds(scores = list(a = s),
                                     split_scores = list(a = cbind(s_h1, s_h2)))
  expect_s3_class(merged, "lpmec_reliability_bounds")
  expect_identical(unname(merged$measure_source["a"]), "scores+split_scores")
  expect_true(is.finite(
    merged$reliability$split_correlation[
      merged$reliability$measure == "a" &
      merged$reliability$design == "pooled"]
  ))
  # Duplicates WITHIN a source remain an error.
  expect_error(
    lpmec_reliability_bounds(scores = list(a = s, a = s)),
    "[Dd]uplicate measure name"
  )
  expect_error(
    lpmec_reliability_bounds(scores = list(a = s), designs = "within"),
    "unit"
  )
  expect_error(
    lpmec_reliability_bounds(scores = list(a = s),
                             unit = rep(1:12, each = 10),
                             designs = "difference"),
    "time"
  )
  expect_error(
    lpmec_reliability_bounds(scores = list(a = s, b = s[-1])),
    "same number of rows"
  )
  expect_error(
    lpmec_reliability_bounds(split_scores = list(a = cbind(s, s, s))),
    "exactly 2 columns"
  )
})
