# Skip all tests on CRAN to avoid timeouts
skip_on_cran()

# Smoke tests for the S3 methods of the 1.2.0 classes: every print/summary
# method must run without error and return invisibly; plots are drawn on a
# null device. Fixtures reuse helper-data.R with tiny n_boot.

expect_prints_invisibly <- function(object) {
  printed <- capture.output(result <- withVisible(print(object)))
  expect_type(printed, "character")
  expect_gt(length(printed), 0L)
  expect_false(result$visible)
  expect_identical(result$value, object)
}

expect_summarizes_invisibly <- function(object) {
  printed <- capture.output(result <- withVisible(summary(object)))
  expect_type(printed, "character")
  expect_gt(length(printed), 0L)
  expect_false(result$visible)
  expect_s3_class(result$value, "data.frame")
  result$value
}

panel_dat <- make_panel_test_data(n_units = 40L, n_periods = 10L,
                                  sig2U = 0.25, seed = 42L)
panel_onerun <- lpmec_panel_onerun(
  Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
  split_scores = list(m1 = cbind(panel_dat$t1, panel_dat$t2)),
  design = "within"
)
panel_agg <- suppressMessages(lpmec_panel(
  Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
  split_scores = list(m1 = cbind(panel_dat$t1, panel_dat$t2)),
  design = "within", n_boot = 2L, seed = 1
))
# scores-only two-measure run: exercises the no-splits fallback paths
panel_scores_agg <- suppressMessages(lpmec_panel(
  Y = panel_dat$Y, unit = panel_dat$unit, time = panel_dat$time,
  scores = list(s1 = panel_dat$t1, s2 = panel_dat$t2),
  design = "within", n_boot = 0L, seed = 2
))

moderator_dat <- make_moderator_test_data(n = 300L, seed = 202L)
moderator_onerun <- suppressMessages(lpmec_moderator_onerun(
  Y = moderator_dat$Y, treatment = moderator_dat$treatment,
  observables = moderator_dat$items, estimation_method = "averaging"
))
moderator_agg <- suppressMessages(lpmec_moderator(
  Y = moderator_dat$Y, treatment = moderator_dat$treatment,
  observables = moderator_dat$items, estimation_method = "averaging",
  n_boot = 2L, n_partition = 1L, seed = 7
))

set.seed(100)
n_bounds <- 200L
bounds_latent <- rnorm(n_bounds)
bounds <- lpmec_reliability_bounds(
  split_scores = list(
    m1 = cbind(bounds_latent + rnorm(n_bounds, sd = 0.6),
               bounds_latent + rnorm(n_bounds, sd = 0.6))
  ),
  scores = list(m2 = bounds_latent + rnorm(n_bounds, sd = 0.7),
                m3 = bounds_latent + rnorm(n_bounds, sd = 0.7))
)

test_that("lpmec_panel_onerun print and summary methods work", {
  expect_prints_invisibly(panel_onerun)
  sum_df <- expect_summarizes_invisibly(panel_onerun)
  expect_equal(row.names(sum_df), panel_onerun$measure_names)
  expect_true(all(c("OLS", "Corrected_OLS", "Corrected_OLS_Lower",
                    "Corrected_OLS_Upper", "Corrected_IV",
                    "Split_Correlation", "Triad", "First_Stage_F") %in%
                    colnames(sum_df)))
})

test_that("lpmec_panel print and summary methods work", {
  expect_prints_invisibly(panel_agg)
  sum_df <- expect_summarizes_invisibly(panel_agg)
  expect_equal(row.names(sum_df), panel_agg$measure_names)
  expect_true(all(c("OLS", "OLS_SE", "Corrected_OLS", "Corrected_OLS_SE",
                    "Corrected_OLS_Lower", "Corrected_OLS_Upper",
                    "Corrected_IV", "Corrected_IV_SE", "Split_Correlation",
                    "Triad", "First_Stage_F") %in% colnames(sum_df)))
  expect_prints_invisibly(panel_scores_agg)
  expect_summarizes_invisibly(panel_scores_agg)
})

test_that("plot.lpmec_panel draws half scores and falls back to score pairs", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(panel_agg))
  expect_silent(plot(panel_scores_agg))
})

test_that("lpmec_moderator_onerun print and summary methods work", {
  expect_prints_invisibly(moderator_onerun)
  sum_df <- expect_summarizes_invisibly(moderator_onerun)
  expect_equal(row.names(sum_df),
               c("Treatment", "Moderator", "Interaction",
                 "Corrected Moderator", "Corrected Interaction"))
  expect_equal(sum_df["Interaction", "Estimate"],
               moderator_onerun$interaction_coef)
})

test_that("lpmec_moderator print and summary methods work", {
  expect_prints_invisibly(moderator_agg)
  sum_df <- expect_summarizes_invisibly(moderator_agg)
  expect_equal(row.names(sum_df),
               c("Treatment", "Moderator", "Interaction",
                 "Corrected Moderator", "Corrected Interaction",
                 "rho_half", "rho_score"))
  expect_true(all(c("Estimate", "SE", "CI_Lower", "CI_Upper") %in%
                    colnames(sum_df)))
})

test_that("plot.lpmec_moderator draws the half-score scatter", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(moderator_agg))
})

test_that("print.lpmec_reliability_bounds works", {
  expect_prints_invisibly(bounds)
  printed <- capture.output(print(bounds))
  expect_true(any(grepl("rho_lo", printed)))
})
