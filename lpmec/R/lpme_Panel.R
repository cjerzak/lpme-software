# lpme_Panel.R -- panel/fixed-effects measurement-error corrections:
# exported lpmec_panel_onerun() and lpmec_panel(), plus the internal
# correction engine .lpmec_panel_corrections(). Ports the verified V2
# reference implementations (estimate_core / run_regressions / corrections)
# onto the shared helpers in lpme_PanelTransforms.R and lpme_Measures.R.

#' Panel corrections from naive fits and reliability estimates
#'
#' Applies the Proposition-3 correction algebra to per-measure naive OLS
#' coefficients and split-IV coefficients estimated under a panel design
#' transform. The OLS correction is
#' \code{b * sqrt(rho_pooled) / rho_design} with split-, triad-, and (for a
#' two-measure system without splits) parallel-pair-based reliability
#' variants. The split-based variant uses the Spearman-Brown stepped-up
#' reliabilities \code{split_rho_pooled}/\code{split_rho_design} so the
#' correction matches the full score entered in the regression (Assumption
#' 3: the splits' error covariance must equal that of the regressed score;
#' raw half-split correlations estimate half-score reliability and would
#' overcorrect); the triad and pair variants are already on the full-score
#' scale. Reliabilities below \code{min_reliability} make the corrected
#' coefficient \code{NA} and are collected into a single warning. The
#' within-measure IV correction multiplies by \code{sqrt} of the raw pooled
#' split correlation (its regressor is a half score) and is never floored
#' (multiplication cannot blow up). The cross-measure IV correction is
#' target-oriented (Proposition 3c): with three or more measures the
#' \code{m_by_l} coefficient is multiplied by \code{sqrt} of target
#' \code{m}'s pooled triad reliability -- never by \code{sqrt} of the
#' pairwise correlation, which is consistent only under equal
#' reliabilities; with exactly two measures the pooled pair correlation is
#' the documented parallel-measures fallback. Design-local corrected
#' coefficients (Proposition 2a) divide \code{ols_coef_local} by \code{sqrt}
#' of the design-scale reliability candidates. When \code{rho_Y} is finite
#' (latent outcome, Proposition 4), every corrected coefficient is
#' additionally divided by \code{sqrt(rho_Y)}.
#'
#' @param ols_coef Named numeric vector of naive design-transformed OLS
#'   coefficients, one per measure.
#' @param iv_coef_a,iv_coef_b Named numeric vectors of within-measure
#'   split-IV coefficients (half 1 instrumented by half 2, and the reverse).
#' @param cross_iv_coef Named numeric vector of cross-measure IV
#'   coefficients with names \code{"<regressor>_by_<instrument>"}.
#' @param measure_names Character vector of measure names.
#' @param split_pooled,split_design Named per-measure raw half-split
#'   correlations under the pooled and estimation designs (used for the IV
#'   multiplication and availability checks).
#' @param split_rho_pooled,split_rho_design Named per-measure Spearman-Brown
#'   stepped-up split reliabilities on the full-score scale (default: the
#'   step-up of the raw correlations); used for the OLS division.
#' @param triad_pooled,triad_design Named per-measure triad reliabilities
#'   under the pooled and estimation designs.
#' @param pair_pooled,pair_design Named per-measure parallel-pair
#'   correlations (cross-measure correlation when exactly two measures are
#'   supplied; \code{NA} otherwise).
#' @param cross_pooled Pooled cross-measure correlation matrix.
#' @param ols_coef_local Optional named numeric vector of naive design-local
#'   OLS coefficients (\code{ols_coef * sd(design-transformed score)}).
#' @param rho_Y Scalar full-score reliability of a latent outcome
#'   (\code{NA} for an observed outcome).
#' @param min_reliability Reliability floor for the OLS division.
#' @param warn Logical; emit the single floor warning.
#'
#' @return List of corrected coefficient vectors and their sources.
#'
#' @noRd
.lpmec_panel_corrections <- function(ols_coef,
                                     iv_coef_a,
                                     iv_coef_b,
                                     cross_iv_coef,
                                     measure_names,
                                     split_pooled,
                                     split_design,
                                     triad_pooled,
                                     triad_design,
                                     pair_pooled,
                                     pair_design,
                                     cross_pooled,
                                     split_rho_pooled = stats::setNames(
                                       .lpmec_spearman_brown(split_pooled),
                                       names(split_pooled)),
                                     split_rho_design = stats::setNames(
                                       .lpmec_spearman_brown(split_design),
                                       names(split_design)),
                                     ols_coef_local = NULL,
                                     rho_Y = NA_real_,
                                     min_reliability,
                                     warn = TRUE) {
  n_measures <- length(measure_names)
  na_measure_vec <- stats::setNames(rep(NA_real_, n_measures), measure_names)
  floored_labels <- character(0L)

  correct_ols <- function(b, rho_pooled, rho_design, label) {
    if (!is.finite(b) || !is.finite(rho_pooled) || !is.finite(rho_design)) {
      return(NA_real_)
    }
    if (rho_pooled < min_reliability || rho_design < min_reliability) {
      floored_labels <<- c(floored_labels, label)
      return(NA_real_)
    }
    b * sqrt(rho_pooled) / rho_design
  }
  iv_factor <- function(rho_pooled) {
    if (is.finite(rho_pooled) && rho_pooled > 0) {
      return(sqrt(rho_pooled))
    }
    NA_real_
  }
  correct_local <- function(b, rho_design, label) {
    if (!is.finite(b) || !is.finite(rho_design)) {
      return(NA_real_)
    }
    if (rho_design < min_reliability) {
      floored_labels <<- c(floored_labels, label)
      return(NA_real_)
    }
    b / sqrt(rho_design)
  }
  # Latent-outcome factor (Prop 4): a further division by sqrt(rho_Y),
  # applied to every corrected estimator at the end. NA rho_Y (observed
  # outcome) leaves the factor at 1; a floored rho_Y makes all corrected
  # coefficients NA.
  outcome_factor <- 1
  if (is.finite(rho_Y)) {
    if (rho_Y < min_reliability) {
      floored_labels <- c(floored_labels, "outcome (rho_Y)")
      outcome_factor <- NA_real_
    } else {
      outcome_factor <- 1 / sqrt(rho_Y)
    }
  }

  corrected_ols_coef_split <- na_measure_vec
  corrected_ols_coef_triad <- na_measure_vec
  corrected_ols_coef_pair <- na_measure_vec
  corrected_ols_coef <- na_measure_vec
  corrected_ols_lower <- na_measure_vec
  corrected_ols_upper <- na_measure_vec
  corrected_ols_source <- stats::setNames(
    rep(NA_character_, n_measures), measure_names
  )

  for (m in measure_names) {
    corrected_ols_coef_split[[m]] <- correct_ols(
      ols_coef[[m]], split_rho_pooled[[m]], split_rho_design[[m]],
      paste0(m, " (split)")
    )
    corrected_ols_coef_triad[[m]] <- correct_ols(
      ols_coef[[m]], triad_pooled[[m]], triad_design[[m]],
      paste0(m, " (triad)")
    )
    # Parallel-pair fallback: with exactly two measures and no within-measure
    # split available, the cross-measure correlation identifies the (common)
    # reliability under a parallel-measures assumption.
    if (n_measures == 2L &&
        !is.finite(split_pooled[[m]]) && !is.finite(split_design[[m]])) {
      corrected_ols_coef_pair[[m]] <- correct_ols(
        ols_coef[[m]], pair_pooled[[m]], pair_design[[m]],
        paste0(m, " (pair)")
      )
    }
    candidates <- c(
      split = corrected_ols_coef_split[[m]],
      triad = corrected_ols_coef_triad[[m]],
      pair = corrected_ols_coef_pair[[m]]
    )
    finite_candidates <- candidates[is.finite(candidates)]
    if (length(finite_candidates) > 0L) {
      corrected_ols_coef[[m]] <- finite_candidates[[1L]]
      corrected_ols_source[[m]] <- names(finite_candidates)[1L]
      corrected_ols_lower[[m]] <- min(finite_candidates)
      corrected_ols_upper[[m]] <- max(finite_candidates)
    }
  }

  # Design-local scale (Prop 2a): the naive local slope (per SD of the
  # transformed score) divided by sqrt of the design-scale reliability.
  if (is.null(ols_coef_local)) {
    ols_coef_local <- na_measure_vec
  }
  corrected_ols_coef_local_split <- na_measure_vec
  corrected_ols_coef_local_triad <- na_measure_vec
  corrected_ols_coef_local_pair <- na_measure_vec
  corrected_ols_coef_local <- na_measure_vec
  corrected_ols_local_lower <- na_measure_vec
  corrected_ols_local_upper <- na_measure_vec
  corrected_ols_local_source <- stats::setNames(
    rep(NA_character_, n_measures), measure_names
  )
  for (m in measure_names) {
    corrected_ols_coef_local_split[[m]] <- correct_local(
      ols_coef_local[[m]], split_rho_design[[m]],
      paste0(m, " (local split)")
    )
    corrected_ols_coef_local_triad[[m]] <- correct_local(
      ols_coef_local[[m]], triad_design[[m]],
      paste0(m, " (local triad)")
    )
    if (n_measures == 2L &&
        !is.finite(split_pooled[[m]]) && !is.finite(split_design[[m]])) {
      corrected_ols_coef_local_pair[[m]] <- correct_local(
        ols_coef_local[[m]], pair_design[[m]],
        paste0(m, " (local pair)")
      )
    }
    local_candidates <- c(
      split = corrected_ols_coef_local_split[[m]],
      triad = corrected_ols_coef_local_triad[[m]],
      pair = corrected_ols_coef_local_pair[[m]]
    )
    finite_local <- local_candidates[is.finite(local_candidates)]
    if (length(finite_local) > 0L) {
      corrected_ols_coef_local[[m]] <- finite_local[[1L]]
      corrected_ols_local_source[[m]] <- names(finite_local)[1L]
      corrected_ols_local_lower[[m]] <- min(finite_local)
      corrected_ols_local_upper[[m]] <- max(finite_local)
    }
  }

  split_iv_factor <- vapply(measure_names, function(m) {
    iv_factor(split_pooled[[m]])
  }, numeric(1L))
  corrected_iv_coef_a <- iv_coef_a * split_iv_factor
  corrected_iv_coef_b <- iv_coef_b * split_iv_factor
  corrected_iv_coef_within <- (corrected_iv_coef_a + corrected_iv_coef_b) / 2

  corrected_cross_iv_coef <- cross_iv_coef
  corrected_cross_iv_pair <- numeric(0L)
  corrected_iv_coef_cross <- na_measure_vec
  if (n_measures >= 2L) {
    # Target-oriented correction (Prop 3c): the m_by_l coefficient targets
    # measure m, so it is multiplied by sqrt of m's pooled triad
    # reliability when 3+ measures identify it. With exactly two measures
    # the pooled pair correlation is the equal-reliability parallel-pair
    # fallback. The pairwise sqrt(r_ml) is never used at M >= 3.
    cross_factor <- function(target, instrument) {
      if (n_measures >= 3L) {
        return(iv_factor(triad_pooled[[target]]))
      }
      iv_factor(cross_pooled[target, instrument])
    }
    for (l in seq_len(n_measures - 1L)) {
      for (m in seq.int(l + 1L, n_measures)) {
        name_l <- measure_names[l]
        name_m <- measure_names[m]
        key_lm <- paste0(name_l, "_by_", name_m)
        key_ml <- paste0(name_m, "_by_", name_l)
        corrected_cross_iv_coef[[key_lm]] <-
          cross_iv_coef[[key_lm]] * cross_factor(name_l, name_m)
        corrected_cross_iv_coef[[key_ml]] <-
          cross_iv_coef[[key_ml]] * cross_factor(name_m, name_l)
        corrected_cross_iv_pair[[paste0(name_l, "_x_", name_m)]] <-
          (corrected_cross_iv_coef[[key_lm]] +
             corrected_cross_iv_coef[[key_ml]]) / 2
      }
    }
    for (m in measure_names) {
      other_measures <- setdiff(measure_names, m)
      corrected_iv_coef_cross[[m]] <- mean(vapply(other_measures, function(l) {
        corrected_cross_iv_coef[[paste0(m, "_by_", l)]]
      }, numeric(1L)))
    }
  }

  corrected_iv_coef <- na_measure_vec
  corrected_iv_source <- stats::setNames(
    rep(NA_character_, n_measures), measure_names
  )
  for (m in measure_names) {
    if (is.finite(corrected_iv_coef_within[[m]])) {
      corrected_iv_coef[[m]] <- corrected_iv_coef_within[[m]]
      corrected_iv_source[[m]] <- "within_split"
    } else if (is.finite(corrected_iv_coef_cross[[m]])) {
      corrected_iv_coef[[m]] <- corrected_iv_coef_cross[[m]]
      corrected_iv_source[[m]] <- "cross_measure"
    }
  }

  # Latent-outcome scaling (Prop 4): a scalar factor common to every
  # corrected estimator (naive coefficients stay naive). NA (floored rho_Y)
  # propagates to every corrected coefficient.
  corrected_ols_coef_split <- corrected_ols_coef_split * outcome_factor
  corrected_ols_coef_triad <- corrected_ols_coef_triad * outcome_factor
  corrected_ols_coef_pair <- corrected_ols_coef_pair * outcome_factor
  corrected_ols_coef <- corrected_ols_coef * outcome_factor
  corrected_ols_lower <- corrected_ols_lower * outcome_factor
  corrected_ols_upper <- corrected_ols_upper * outcome_factor
  corrected_ols_coef_local_split <-
    corrected_ols_coef_local_split * outcome_factor
  corrected_ols_coef_local_triad <-
    corrected_ols_coef_local_triad * outcome_factor
  corrected_ols_coef_local_pair <-
    corrected_ols_coef_local_pair * outcome_factor
  corrected_ols_coef_local <- corrected_ols_coef_local * outcome_factor
  corrected_ols_local_lower <- corrected_ols_local_lower * outcome_factor
  corrected_ols_local_upper <- corrected_ols_local_upper * outcome_factor
  corrected_iv_coef_a <- corrected_iv_coef_a * outcome_factor
  corrected_iv_coef_b <- corrected_iv_coef_b * outcome_factor
  corrected_iv_coef_within <- corrected_iv_coef_within * outcome_factor
  corrected_cross_iv_coef <- corrected_cross_iv_coef * outcome_factor
  corrected_cross_iv_pair <- corrected_cross_iv_pair * outcome_factor
  corrected_iv_coef_cross <- corrected_iv_coef_cross * outcome_factor
  corrected_iv_coef <- corrected_iv_coef * outcome_factor

  if (warn && length(floored_labels) > 0L) {
    warning(
      "Reliability estimate(s) below 'min_reliability' (", min_reliability,
      "); corrected coefficient(s) reported as NA for: ",
      paste(unique(floored_labels), collapse = ", "), ".",
      call. = FALSE
    )
  }

  list(
    corrected_ols_coef = corrected_ols_coef,
    corrected_ols_coef_split = corrected_ols_coef_split,
    corrected_ols_coef_triad = corrected_ols_coef_triad,
    corrected_ols_coef_pair = corrected_ols_coef_pair,
    corrected_ols_lower = corrected_ols_lower,
    corrected_ols_upper = corrected_ols_upper,
    corrected_ols_source = corrected_ols_source,
    corrected_ols_coef_local = corrected_ols_coef_local,
    corrected_ols_coef_local_split = corrected_ols_coef_local_split,
    corrected_ols_coef_local_triad = corrected_ols_coef_local_triad,
    corrected_ols_coef_local_pair = corrected_ols_coef_local_pair,
    corrected_ols_local_lower = corrected_ols_local_lower,
    corrected_ols_local_upper = corrected_ols_local_upper,
    corrected_ols_local_source = corrected_ols_local_source,
    corrected_iv_coef_a = corrected_iv_coef_a,
    corrected_iv_coef_b = corrected_iv_coef_b,
    corrected_iv_coef_within = corrected_iv_coef_within,
    corrected_cross_iv_coef = corrected_cross_iv_coef,
    corrected_cross_iv_pair = corrected_cross_iv_pair,
    corrected_iv_coef_cross = corrected_iv_coef_cross,
    corrected_iv_coef = corrected_iv_coef,
    corrected_iv_source = corrected_iv_source
  )
}

#' Single-run panel/fixed-effects measurement-error correction
#'
#' Estimates the effect of a latent predictor on a panel outcome under a
#' fixed-effects or differencing design, and corrects the naive coefficient
#' for the design-amplified attenuation caused by measurement error in the
#' latent-variable scores (Propositions 1, 3, and 7 of the accompanying
#' working paper; Proposition 2a for the design-local scale and
#' Proposition 4 for latent outcomes). Latent measures may be supplied as
#' raw item batteries (\code{observables}, scored via
#' \code{\link{lpmec_onerun}}), as pre-computed half scores
#' (\code{split_scores}), or as full pre-computed scores (\code{scores});
#' the three sources are merged by measure name.
#'
#' @param Y Numeric outcome vector (one value per panel row). Non-finite
#'   values are treated as missing. Optional when \code{Y_split_scores} is
#'   supplied (the outcome score is then built from the two halves).
#' @param unit Vector of unit (cluster) identifiers, one per row. Required
#'   for the \code{"within"}, \code{"twoway"}, and \code{"difference"}
#'   designs (it defines the fixed-effect groups, the exact-gap
#'   differencing paths, and the cluster-robust standard errors). May be
#'   \code{NULL} for \code{design = "pooled"} (each row is then its own
#'   cluster, covering the pure cross-sectional case).
#' @param time Optional numeric vector of time identifiers, one per row.
#'   Required for the \code{"twoway"} and \code{"difference"} designs.
#'   (unit, time) pairs must be unique.
#' @param observables Optional list of item matrices or data frames, one per
#'   measure (a single matrix is treated as one measure). Each measure is
#'   scored via \code{\link{lpmec_onerun}} with \code{estimation_method} and
#'   must have at least 4 columns; measures with fewer components must be
#'   supplied through \code{split_scores}.
#' @param scores Optional named list (or matrix/data frame with one column
#'   per measure) of pre-computed full measure scores.
#' @param split_scores Optional named list of n x 2 matrices of half scores,
#'   one per measure.
#' @param covariates Optional matrix or data frame of observed covariates;
#'   they are design-transformed together with the outcome and scores and
#'   enter both the OLS and both IV stages. Values must be finite.
#' @param Y_split_scores Optional n x 2 matrix of outcome half scores for a
#'   \emph{latent} outcome (Proposition 4). Each half is pooled z-scored;
#'   the outcome reliability \code{rho_Y} is the Spearman-Brown step-up of
#'   the half correlation, and every corrected estimator is additionally
#'   divided by \code{sqrt(rho_Y)}. When \code{Y} is not supplied, the
#'   outcome score is the pooled z-score of the mean of the two z-scored
#'   halves; when \code{Y} is also supplied, it is used as the outcome
#'   score and the halves only estimate \code{rho_Y} (a sensitivity
#'   calculation unless the halves' error covariance matches that of
#'   \code{Y}).
#' @param design Panel design transform for estimation: one of
#'   \code{"within"} (unit demeaning), \code{"twoway"} (iterated
#'   unit-and-time demeaning), \code{"difference"} (exact-gap k-period
#'   differencing), or \code{"pooled"} (identity). Default \code{"within"}.
#' @param diff_k Positive integer gap for the \code{"difference"} design.
#'   Rows lacking a same-unit observation at exactly \code{time - diff_k}
#'   are set to \code{NA} rather than differenced against a nearer one.
#' @param estimation_method Estimation method(s) passed to
#'   \code{\link{lpmec_onerun}} for \code{observables}-mode measures (single
#'   value or one per measure). Default \code{"averaging"}.
#' @param min_reliability Reliability floor for the corrected-OLS division.
#'   Estimated reliabilities below this value (or negative) yield \code{NA}
#'   corrected OLS coefficients with a single warning. Default \code{0.05}.
#' @param min_cor_n Minimum number of complete pairs required to report any
#'   correlation. Default \code{30}.
#' @param demean_iterations Number of iterated demeaning passes for the
#'   \code{"twoway"} design (one pass is exact only on balanced panels).
#'   Default \code{25}.
#' @param partition Optional named list of fixed split-half partitions,
#'   keyed by measure name, passed to \code{\link{lpmec_onerun}} for
#'   \code{observables}-mode measures.
#' @param ... Additional arguments passed to \code{\link{lpmec_onerun}} for
#'   \code{observables}-mode measures (e.g., \code{ordinal},
#'   \code{mcmc_control}).
#'
#' @return A list of class \code{lpmec_panel_onerun} containing, with one
#' named entry per measure unless noted:
#' \itemize{
#'   \item \code{ols_coef}, \code{ols_se}, \code{ols_n_obs},
#'     \code{ols_n_clusters}: naive design-transformed OLS coefficient of
#'     \code{Y} on the measure score with cluster-robust (by unit) standard
#'     error.
#'   \item \code{corrected_ols_coef}, \code{corrected_ols_coef_split},
#'     \code{corrected_ols_coef_triad}, \code{corrected_ols_coef_pair},
#'     \code{corrected_ols_lower}, \code{corrected_ols_upper},
#'     \code{corrected_ols_source}: corrected OLS coefficient(s)
#'     \code{b * sqrt(rho_pooled) / rho_design} per reliability variant
#'     (the split variant uses the Spearman-Brown stepped-up reliabilities
#'     \code{split_rho_score}/\code{design_split_rho_score}), the headline
#'     value (split when available, else triad, else the two-measure
#'     parallel pair), and the [min, max] sensitivity range over the finite
#'     variants.
#'   \item \code{sd_design_x}, \code{ols_coef_local},
#'     \code{corrected_ols_coef_local(_split/_triad/_pair)},
#'     \code{corrected_ols_local_lower/_upper},
#'     \code{corrected_ols_local_source}: design-local-scale quantities
#'     (Proposition 2a): the sd of the design-transformed score on the
#'     estimation sample, the naive local slope
#'     \code{ols_coef * sd_design_x} (effect per SD of the
#'     \emph{transformed} latent trait), and its corrections
#'     \code{ols_coef_local / sqrt(rho_design)} per design-scale
#'     reliability variant. The local estimand moves with the design and is
#'     not comparable across specifications.
#'   \item \code{split_correlation}, \code{design_split_correlation},
#'     \code{split_rho_score}, \code{design_split_rho_score}, \code{triad},
#'     \code{design_triad}, \code{pair_correlation},
#'     \code{design_pair_correlation}: reliability estimates under the
#'     pooled and estimation designs; \code{*_rho_score} are the
#'     Spearman-Brown step-ups of the raw half-split correlations to the
#'     full-score scale.
#'   \item \code{rho_Y_half}, \code{rho_Y}: with \code{Y_split_scores}, the
#'     raw outcome half-score correlation and its Spearman-Brown step-up
#'     (\code{NA} for an observed outcome).
#'   \item \code{iv_coef_a}, \code{iv_coef_b}, \code{iv_coef},
#'     \code{iv_se_a}, \code{iv_se_b}, \code{first_stage_fstat_a},
#'     \code{first_stage_fstat_b}: within-measure split-IV fits (half 1
#'     instrumented by half 2 and the reverse) on the design-transformed
#'     data, with clustered first-stage F statistics.
#'   \item \code{corrected_iv_coef_a}, \code{corrected_iv_coef_b},
#'     \code{corrected_iv_coef_within}: within-measure IV corrected by
#'     \code{sqrt} of the pooled split correlation.
#'   \item \code{cross_iv_coef}, \code{cross_iv_se},
#'     \code{cross_first_stage_fstat}, \code{corrected_cross_iv_coef}
#'     (named \code{"<regressor>_by_<instrument>"}),
#'     \code{corrected_cross_iv_pair} (direction-averaged, named
#'     \code{"<a>_x_<b>"}), \code{corrected_iv_coef_cross} (per regressor,
#'     averaged over instruments): cross-measure IV fits with the
#'     target-oriented correction of Proposition 3c -- \code{sqrt} of the
#'     regressor (target) measure's pooled triad reliability when 3+
#'     measures are supplied, or \code{sqrt} of the pooled pair correlation
#'     under the two-measure parallel-measures fallback.
#'   \item \code{corrected_iv_coef}, \code{corrected_iv_source},
#'     \code{first_stage_fstat}: headline corrected IV per measure
#'     (within-measure when available, else cross-measure) and the minimum
#'     first-stage F statistic among the directions it averages.
#'   \item \code{reliability}: data frame with one row per (design, measure)
#'     holding \code{split_correlation}, \code{split_n}, \code{rho_split},
#'     \code{triad}, \code{rho_lo}, \code{rho_hi} (the Proposition-7
#'     sensitivity range).
#'   \item \code{cor_matrices}, \code{cor_n_matrices}: per-design
#'     cross-measure correlation matrices and complete-pair counts.
#'   \item \code{x_est}, \code{x_est1}, \code{x_est2}: pooled z-scored,
#'     sign-aligned measure scores and half scores; \code{scalar_runs}: the
#'     underlying \code{\link{lpmec_onerun}} fits for
#'     \code{observables}-mode measures.
#'   \item Metadata: \code{measure_names}, \code{n_measures},
#'     \code{measure_source}, \code{has_splits}, \code{sign_flipped},
#'     \code{covariate_names}, \code{design}, \code{diff_k}, \code{n_obs},
#'     \code{n_units}, \code{min_reliability}, \code{min_cor_n},
#'     \code{demean_iterations}.
#' }
#'
#' @details
#' Under the panel design transform, classical measurement error attenuates
#' the naive coefficient by the design reliability \code{rho_design} of the
#' (pooled-standardized) score, while the pooled standardization inflates it
#' by \code{1 / sqrt(rho_pooled)}, so
#' \code{plim b = beta * sqrt(rho_pooled) * rho_design / rho_pooled}
#' simplifies to \code{beta * rho_design / sqrt(rho_pooled)} and the
#' corrected OLS coefficient is \code{b * sqrt(rho_pooled) / rho_design}.
#' Fixed-effects and differencing transforms strip the persistent part of
#' the latent signal but not the transient measurement error, so
#' \code{rho_design} is typically far below \code{rho_pooled} and pooled
#' corrections badly understate the bias.
#'
#' \strong{Score-scale matching (Assumption 3).} The reliabilities entering
#' the OLS division must be those of the score actually regressed. The raw
#' half-split correlation estimates the reliability of a \emph{half} score,
#' while the regression uses the full score (the average of the two z-scored
#' halves, or the full-battery estimate), so the split-based correction uses
#' the Spearman-Brown step-up \code{2r / (1 + r)} of both the pooled and the
#' design split correlations. The step-up is exact when the full score is
#' the average of two independent parallel halves and a parallel-forms
#' approximation for full-battery estimates. When a measure merges a
#' user-supplied full score with separate half scores
#' (\code{scores} + \code{split_scores}), the stepped-up value matches the
#' supplied score's error covariance only under that same parallel-forms
#' assumption, so the split-based correction is then a sensitivity
#' calculation.
#'
#' \strong{IV multiplies, OLS divides.} The corrected OLS coefficient
#' divides by the design reliability, which is estimated with noise and can
#' be near zero under aggressive transforms; corrections with
#' \code{rho_design < min_reliability} (or with the pooled reliability below
#' the floor) are therefore reported as \code{NA} rather than divided out.
#' The split-IV correction instead multiplies the design-transformed IV
#' coefficient by \code{sqrt(rho_pooled)}: the instrument (the other half
#' score or another measure) removes the design-reliability denominator, so
#' the corrected IV estimate stays stable exactly where the OLS correction
#' becomes fragile. Do not swap the two operations.
#'
#' Reliability variants for the OLS correction: the Spearman-Brown
#' stepped-up within-measure split reliability (exceeds construct-relevant
#' reliability when the halves share systematic error, Proposition 7a), the
#' cross-measure triad \code{r_ml * r_mk / r_lk} (requires 3+ measures;
#' consistent under pairwise-orthogonal total errors, Proposition 7b, and a
#' lower bound only under the target-specific ratio condition of
#' Proposition 7c -- shared error can move it in either direction), and --
#' only for a two-measure system where a measure has no split halves -- the
#' cross-measure pair correlation, which identifies the common reliability
#' under a parallel-measures assumption. \code{corrected_ols_lower} and
#' \code{corrected_ols_upper} span the finite variants; this is a
#' sensitivity range, not a partial-identification interval, unless the
#' corresponding Proposition 7 condition is maintained.
#'
#' \strong{Covariate caveat.} With covariates, the scalar OLS correction is
#' exact only when the covariates are uncorrelated with the latent
#' predictor; otherwise the measurement error also contaminates the
#' covariate coefficients and the correction of the latent slope is
#' approximate. The split-IV estimator remains consistent for the latent
#' slope in the presence of exogenous covariates, which enter both stages.
#'
#' @examples
#' \donttest{
#' set.seed(100)
#' n_units <- 60
#' n_periods <- 12
#' unit <- rep(seq_len(n_units), times = n_periods)
#' time <- rep(seq_len(n_periods), each = n_units)
#' latent <- rnorm(n_units)[unit] + rnorm(n_units * n_periods, sd = 0.7)
#' half1 <- latent + rnorm(n_units * n_periods, sd = 0.5)
#' half2 <- latent + rnorm(n_units * n_periods, sd = 0.5)
#' Y <- 0.4 * latent + rnorm(n_units)[unit] +
#'   rnorm(n_units * n_periods)
#'
#' run <- lpmec_panel_onerun(
#'   Y = Y, unit = unit, time = time,
#'   split_scores = list(m1 = cbind(half1, half2)),
#'   design = "within"
#' )
#' run$ols_coef
#' run$corrected_ols_coef
#' run$corrected_iv_coef
#' }
#'
#' @importFrom sandwich vcovCL
#'
#' @export
lpmec_panel_onerun <- function(Y = NULL,
                               unit = NULL,
                               time = NULL,
                               observables = NULL,
                               scores = NULL,
                               split_scores = NULL,
                               covariates = NULL,
                               Y_split_scores = NULL,
                               design = c("within", "twoway", "difference", "pooled"),
                               diff_k = 1L,
                               estimation_method = "averaging",
                               min_reliability = 0.05,
                               min_cor_n = 30L,
                               demean_iterations = 25L,
                               partition = NULL,
                               ...) {
  design <- match.arg(design)

  # Latent outcome (Prop 4): resolve the outcome score and its reliability
  rho_Y_half <- rho_Y <- NA_real_
  if (!is.null(Y_split_scores)) {
    y_halves <- .lpmec_numeric_observable_matrix(
      Y_split_scores, label = "'Y_split_scores'"
    )
    if (ncol(y_halves) != 2L) {
      stop("'Y_split_scores' must have exactly 2 columns ",
           "(one per outcome half score).")
    }
    y_half1 <- .lpmec_zscore(y_halves[, 1L])
    y_half2 <- .lpmec_zscore(y_halves[, 2L])
    rho_Y_half <- suppressWarnings(
      stats::cor(y_half1, y_half2, use = "pairwise.complete.obs")
    )
    if (!is.finite(rho_Y_half)) {
      stop("Could not estimate the outcome reliability: the correlation ",
           "of the 'Y_split_scores' halves is not finite.")
    }
    rho_Y <- as.numeric(.lpmec_spearman_brown(rho_Y_half))
    if (is.null(Y)) {
      Y <- .lpmec_zscore(rowMeans(cbind(y_half1, y_half2)))
    } else if (length(Y) != nrow(y_halves)) {
      stop("'Y_split_scores' must have one row per element of 'Y'.")
    }
  }
  if (is.null(Y)) {
    stop("'Y' is required unless 'Y_split_scores' is supplied.")
  }
  if (!is.numeric(Y)) {
    stop("'Y' must be a numeric vector.")
  }
  Y <- as.numeric(Y)
  Y[!is.finite(Y)] <- NA_real_

  prep <- .lpmec_prepare_panel_inputs(
    n_obs = length(Y),
    unit = unit,
    time = time,
    design = design,
    diff_k = diff_k,
    min_reliability = min_reliability,
    min_cor_n = min_cor_n,
    demean_iterations = demean_iterations
  )

  measures <- .lpmec_resolve_measures(
    observables = observables,
    scores = scores,
    split_scores = split_scores,
    Y = Y,
    estimation_method = estimation_method,
    partitions = partition,
    ...
  )
  if (measures$n_obs != prep$n_obs) {
    stop("'Y' must be a numeric vector with one value per row of the ",
         "measure inputs. Received ", prep$n_obs, " outcome value(s) for ",
         measures$n_obs, " measure row(s).")
  }
  measure_names <- measures$measure_names
  n_measures <- measures$n_measures
  if (n_measures == 1L && !any(measures$has_splits)) {
    stop("A single measure without split halves is not identified: supply ",
         "the measure via 'observables' or 'split_scores', or provide at ",
         "least two measures.")
  }

  covariate_matrix <- .lpmec_prepare_covariates(covariates, prep$n_obs)

  # Reliability table under the pooled and estimation designs
  rel_designs <- unique(c("pooled", design))
  rel <- .lpmec_reliability_table(
    measures = measures,
    unit = prep$unit,
    time = prep$time,
    designs = rel_designs,
    diff_k = prep$diff_k,
    demean_iterations = prep$demean_iterations,
    min_cor_n = prep$min_cor_n
  )
  reliability <- .lpmec_bounds_from_reliability(
    rel$table,
    min_reliability = prep$min_reliability,
    warn = FALSE
  )

  rel_lookup <- function(design_name, field) {
    rows <- reliability[reliability$design == design_name, , drop = FALSE]
    stats::setNames(
      rows[[field]][match(measure_names, rows$measure)],
      measure_names
    )
  }
  split_correlation <- rel_lookup("pooled", "split_correlation")
  design_split_correlation <- rel_lookup(design, "split_correlation")
  split_rho_score <- rel_lookup("pooled", "rho_split")
  design_split_rho_score <- rel_lookup(design, "rho_split")
  triad <- rel_lookup("pooled", "triad")
  design_triad <- rel_lookup(design, "triad")
  na_measure_vec <- stats::setNames(rep(NA_real_, n_measures), measure_names)
  pair_correlation <- na_measure_vec
  design_pair_correlation <- na_measure_vec
  if (n_measures == 2L) {
    pair_correlation[] <- rel$cor_matrices[["pooled"]][1L, 2L]
    design_pair_correlation[] <- rel$cor_matrices[[design]][1L, 2L]
  }

  # Design transforms of the outcome, scores, half scores, and covariates
  transform_columns <- function(x) {
    .lpmec_panel_transform(
      x,
      unit = prep$unit,
      time = prep$time,
      design = design,
      diff_k = prep$diff_k,
      demean_iterations = prep$demean_iterations
    )
  }
  Y_design <- transform_columns(Y)
  x_design <- transform_columns(measures$x_est)
  x1_design <- transform_columns(measures$x_est1)
  x2_design <- transform_columns(measures$x_est2)
  covariates_design <- if (ncol(covariate_matrix) > 0L) {
    transform_columns(covariate_matrix)
  } else {
    NULL
  }

  # Naive cluster-robust OLS per measure
  ols_coef <- ols_se <- sd_design_x <- na_measure_vec
  ols_n_obs <- ols_n_clusters <- stats::setNames(
    rep(NA_integer_, n_measures), measure_names
  )
  for (m in measure_names) {
    ols_fit <- .lpmec_panel_ols(
      Y_design, x_design[, m],
      covariates = covariates_design, cluster = prep$unit
    )
    ols_coef[[m]] <- ols_fit$coef
    ols_se[[m]] <- ols_fit$se
    ols_n_obs[[m]] <- as.integer(ols_fit$n_obs)
    ols_n_clusters[[m]] <- as.integer(ols_fit$n_clusters)
    sd_design_x[[m]] <- ols_fit$x_sd
  }
  # Design-local scale (Prop 2a): the slope of the design-transformed
  # outcome on the re-standardized transformed score, exact in sample
  # because rescaling a regressor by a scalar rescales its coefficient.
  ols_coef_local <- ols_coef * sd_design_x

  # Within-measure split-IV (both directions) per measure with half scores
  iv_coef_a <- iv_coef_b <- iv_se_a <- iv_se_b <- na_measure_vec
  first_stage_fstat_a <- first_stage_fstat_b <- na_measure_vec
  for (m in measure_names) {
    if (!measures$has_splits[[m]]) {
      next
    }
    iv_fit_a <- .lpmec_panel_iv(
      Y_design, x1_design[, m], x2_design[, m],
      covariates = covariates_design, cluster = prep$unit
    )
    iv_fit_b <- .lpmec_panel_iv(
      Y_design, x2_design[, m], x1_design[, m],
      covariates = covariates_design, cluster = prep$unit
    )
    iv_coef_a[[m]] <- iv_fit_a$coef
    iv_coef_b[[m]] <- iv_fit_b$coef
    iv_se_a[[m]] <- iv_fit_a$se
    iv_se_b[[m]] <- iv_fit_b$se
    first_stage_fstat_a[[m]] <- iv_fit_a$first_stage_fstat
    first_stage_fstat_b[[m]] <- iv_fit_b$first_stage_fstat
  }
  iv_coef <- (iv_coef_a + iv_coef_b) / 2

  # Cross-measure split-IV over all ordered (regressor, instrument) pairs
  cross_names <- character(0L)
  if (n_measures >= 2L) {
    for (a in measure_names) {
      for (b in setdiff(measure_names, a)) {
        cross_names <- c(cross_names, paste0(a, "_by_", b))
      }
    }
  }
  cross_iv_coef <- cross_iv_se <- cross_first_stage_fstat <-
    stats::setNames(rep(NA_real_, length(cross_names)), cross_names)
  if (n_measures >= 2L) {
    for (a in measure_names) {
      for (b in setdiff(measure_names, a)) {
        cross_fit <- .lpmec_panel_iv(
          Y_design, x_design[, a], x_design[, b],
          covariates = covariates_design, cluster = prep$unit
        )
        key <- paste0(a, "_by_", b)
        cross_iv_coef[[key]] <- cross_fit$coef
        cross_iv_se[[key]] <- cross_fit$se
        cross_first_stage_fstat[[key]] <- cross_fit$first_stage_fstat
      }
    }
  }

  corrections <- .lpmec_panel_corrections(
    ols_coef = ols_coef,
    iv_coef_a = iv_coef_a,
    iv_coef_b = iv_coef_b,
    cross_iv_coef = cross_iv_coef,
    measure_names = measure_names,
    split_pooled = split_correlation,
    split_design = design_split_correlation,
    split_rho_pooled = split_rho_score,
    split_rho_design = design_split_rho_score,
    triad_pooled = triad,
    triad_design = design_triad,
    pair_pooled = pair_correlation,
    pair_design = design_pair_correlation,
    cross_pooled = rel$cor_matrices[["pooled"]],
    ols_coef_local = ols_coef_local,
    rho_Y = rho_Y,
    min_reliability = prep$min_reliability,
    warn = TRUE
  )

  # Headline first-stage F: minimum over the directions the headline IV uses
  first_stage_fstat <- na_measure_vec
  for (m in measure_names) {
    source_m <- corrections$corrected_iv_source[[m]]
    if (is.na(source_m)) {
      next
    }
    fstats_m <- if (source_m == "within_split") {
      c(first_stage_fstat_a[[m]], first_stage_fstat_b[[m]])
    } else {
      cross_first_stage_fstat[
        paste0(m, "_by_", setdiff(measure_names, m))
      ]
    }
    fstats_m <- fstats_m[is.finite(fstats_m)]
    if (length(fstats_m) > 0L) {
      first_stage_fstat[[m]] <- min(fstats_m)
    }
  }

  results <- c(
    list(
      measure_names = measure_names,
      n_measures = n_measures,
      measure_source = measures$source,
      has_splits = measures$has_splits,
      sign_flipped = measures$sign_flipped,
      covariate_names = colnames(covariate_matrix),
      design = design,
      diff_k = prep$diff_k,
      n_obs = prep$n_obs,
      n_units = length(unique(prep$unit)),
      ols_coef = ols_coef,
      ols_se = ols_se,
      ols_n_obs = ols_n_obs,
      ols_n_clusters = ols_n_clusters,
      sd_design_x = sd_design_x,
      ols_coef_local = ols_coef_local,
      split_correlation = split_correlation,
      design_split_correlation = design_split_correlation,
      split_rho_score = split_rho_score,
      design_split_rho_score = design_split_rho_score,
      triad = triad,
      design_triad = design_triad,
      pair_correlation = pair_correlation,
      design_pair_correlation = design_pair_correlation,
      rho_Y_half = rho_Y_half,
      rho_Y = rho_Y,
      iv_coef_a = iv_coef_a,
      iv_coef_b = iv_coef_b,
      iv_coef = iv_coef,
      iv_se_a = iv_se_a,
      iv_se_b = iv_se_b,
      first_stage_fstat_a = first_stage_fstat_a,
      first_stage_fstat_b = first_stage_fstat_b,
      cross_iv_coef = cross_iv_coef,
      cross_iv_se = cross_iv_se,
      cross_first_stage_fstat = cross_first_stage_fstat,
      first_stage_fstat = first_stage_fstat
    ),
    corrections,
    list(
      reliability = reliability,
      cor_matrices = rel$cor_matrices,
      cor_n_matrices = rel$cor_n_matrices,
      x_est = measures$x_est,
      x_est1 = measures$x_est1,
      x_est2 = measures$x_est2,
      scalar_runs = measures$scalar_runs,
      min_reliability = prep$min_reliability,
      min_cor_n = prep$min_cor_n,
      demean_iterations = prep$demean_iterations
    )
  )
  class(results) <- "lpmec_panel_onerun"
  results
}

# .lpmec_aggregate_by_boot() requires two or more columns (its vapply/t()
# combination drops the matrix structure for a single column); this wrapper
# delegates to it and handles the single-measure column directly.
.lpmec_panel_aggregate_by_boot <- function(values, boot_ids, aggregation_fn) {
  if (ncol(values) != 1L) {
    return(.lpmec_aggregate_by_boot(values, boot_ids, aggregation_fn))
  }
  boots <- sort(unique(boot_ids))
  matrix(
    vapply(boots, function(boot) {
      aggregation_fn(values[boot_ids == boot, 1L])
    }, numeric(1L)),
    ncol = 1L,
    dimnames = list(NULL, colnames(values))
  )
}

#' Aggregated panel/fixed-effects measurement-error correction
#'
#' Runs \code{\link{lpmec_panel_onerun}} on the original sample and on
#' cluster (unit) bootstrap resamples, optionally over repeated split-half
#' partitions for \code{observables}-mode measures, and aggregates the
#' naive, corrected, and reliability quantities across runs. Units are
#' resampled with replacement and each draw receives a fresh pseudo-id, so a
#' unit drawn twice enters as two distinct clusters/fixed-effect groups; the
#' whole pipeline (scoring, sign alignment, transforms, reliabilities,
#' regressions, corrections) is re-run on every (bootstrap, partition) pair.
#'
#' @inheritParams lpmec_panel_onerun
#' @param n_boot Non-negative integer number of cluster-bootstrap
#'   replications. Default \code{32}.
#' @param n_partition Positive integer number of split-half partitions per
#'   original or bootstrap sample. Partitions only vary for
#'   \code{observables}-mode measures; without \code{observables} the value
#'   is coerced to 1 with a message. Default \code{10}.
#' @param partition_aggregation Aggregation strategy across partitions
#'   within each bootstrap replication: \code{"median"},
#'   \code{"winsorized_mean"}, \code{"trimmed_mean"}, or a function. See
#'   \code{\link{lpmec}}.
#' @param partition_aggregation_probs Quantile probabilities for winsorized
#'   or trimmed partition aggregation.
#' @param return_intermediaries Logical. If \code{TRUE}, returns the
#'   per-run coefficient matrices as \code{Intermediary_*} entries along
#'   with \code{Intermediary_BootIndex} and
#'   \code{Intermediary_PartitionIndex}.
#' @param seed Optional seed applied locally to the scoring partitions and
#'   the bootstrap (the caller's random-number state is restored on exit).
#'
#' @return A list of class \code{lpmec_panel}. For each aggregated quantity
#' \code{<field>} of \code{\link{lpmec_panel_onerun}} -- \code{ols_coef},
#' \code{corrected_ols_coef(_split/_triad/_pair)},
#' \code{corrected_ols_lower/upper}, \code{sd_design_x},
#' \code{ols_coef_local},
#' \code{corrected_ols_coef_local(_split/_triad/_pair)},
#' \code{corrected_ols_local_lower/upper}, \code{iv_coef},
#' \code{corrected_iv_coef(_within/_cross)}, \code{split_correlation},
#' \code{design_split_correlation}, \code{split_rho_score},
#' \code{design_split_rho_score}, \code{triad}, \code{design_triad},
#' \code{pair_correlation}, \code{design_pair_correlation},
#' \code{rho_Y_half}, \code{rho_Y},
#' \code{first_stage_fstat}, \code{cross_iv_coef},
#' \code{corrected_cross_iv_coef}, \code{corrected_cross_iv_pair} -- the
#' result holds the original-sample partition-aggregated estimate
#' \code{<field>} plus bootstrap \code{<field>_se}, \code{<field>_lower},
#' and \code{<field>_upper} (percentile 2.5\%/97.5\%). Also included:
#' \code{ols_cluster_se} (analytic cluster-robust SE from the original
#' sample), the original-sample \code{reliability} table,
#' \code{cor_matrices}, \code{cor_n_matrices}, \code{x_est}, \code{x_est1},
#' \code{x_est2}, \code{corrected_ols_source},
#' \code{corrected_ols_local_source}, \code{corrected_iv_source},
#' metadata (\code{measure_names}, \code{n_measures}, \code{design},
#' \code{diff_k}, \code{n_obs}, \code{n_units}, \code{n_boot},
#' \code{n_partition}, \code{boot_n_failed}), and -- when
#' \code{return_intermediaries = TRUE} -- \code{Intermediary_BootIndex},
#' \code{Intermediary_PartitionIndex}, and one \code{Intermediary_<field>}
#' matrix of per-run values (rows ordered as the runs; bootstrap runs that
#' failed are \code{NA} rows).
#'
#' @details
#' Aggregation follows the package convention: per-run values are first
#' aggregated across partitions within each bootstrap replication
#' (\code{partition_aggregation}), then the original-sample aggregate is
#' the point estimate and the bootstrap aggregates provide the standard
#' error and percentile interval. Bootstrap replications that fail are
#' caught, recorded as \code{NA} rows, counted in \code{boot_n_failed}, and
#' excluded from the bootstrap summaries.
#'
#' See \code{\link{lpmec_panel_onerun}} for the correction algebra,
#' including the "IV multiplies, OLS divides" asymmetry and the covariate
#' caveat.
#'
#' @examples
#' \donttest{
#' set.seed(100)
#' n_units <- 60
#' n_periods <- 12
#' unit <- rep(seq_len(n_units), times = n_periods)
#' time <- rep(seq_len(n_periods), each = n_units)
#' latent <- rnorm(n_units)[unit] + rnorm(n_units * n_periods, sd = 0.7)
#' half1 <- latent + rnorm(n_units * n_periods, sd = 0.5)
#' half2 <- latent + rnorm(n_units * n_periods, sd = 0.5)
#' Y <- 0.4 * latent + rnorm(n_units)[unit] +
#'   rnorm(n_units * n_periods)
#'
#' result <- lpmec_panel(
#'   Y = Y, unit = unit, time = time,
#'   split_scores = list(m1 = cbind(half1, half2)),
#'   design = "within",
#'   n_boot = 4L,
#'   seed = 1
#' )
#' result$corrected_ols_coef
#' result$corrected_ols_coef_se
#' }
#'
#' @export
lpmec_panel <- function(Y = NULL,
                        unit = NULL,
                        time = NULL,
                        observables = NULL,
                        scores = NULL,
                        split_scores = NULL,
                        covariates = NULL,
                        Y_split_scores = NULL,
                        design = c("within", "twoway", "difference", "pooled"),
                        diff_k = 1L,
                        estimation_method = "averaging",
                        min_reliability = 0.05,
                        min_cor_n = 30L,
                        demean_iterations = 25L,
                        n_boot = 32L,
                        n_partition = 10L,
                        partition_aggregation = "median",
                        partition_aggregation_probs = c(0.01, 0.99),
                        return_intermediaries = TRUE,
                        seed = NULL,
                        ...) {
  design <- match.arg(design)
  if (missing(Y)) {
    Y <- NULL
  }
  if (is.null(Y) && is.null(Y_split_scores)) {
    stop("'Y' is required unless 'Y_split_scores' is supplied.")
  }
  if (!is.null(Y) && !is.numeric(Y)) {
    stop("'Y' must be a numeric vector.")
  }
  n_obs_input <- if (!is.null(Y)) length(Y) else nrow(as.matrix(Y_split_scores))
  if (!is.null(Y_split_scores) &&
      nrow(as.matrix(Y_split_scores)) != n_obs_input) {
    stop("'Y_split_scores' must have one row per element of 'Y'.")
  }
  if (missing(unit)) {
    unit <- NULL
  }
  if (is.null(unit) && design != "pooled") {
    stop("'unit' is required for design(s): ", design, ".")
  }
  if (!is.null(unit) && length(unit) != n_obs_input) {
    stop("'unit' must have the same length as 'Y'.")
  }
  if (!is.null(time) && length(time) != n_obs_input) {
    stop("'time' must have the same length as 'Y'.")
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || !is.finite(n_boot) ||
      n_boot != floor(n_boot) || n_boot < 0) {
    stop("'n_boot' must be a single non-negative integer.")
  }
  if (!is.numeric(n_partition) || length(n_partition) != 1L ||
      !is.finite(n_partition) || n_partition != floor(n_partition) ||
      n_partition < 1) {
    stop("'n_partition' must be a single positive integer.")
  }
  n_boot <- as.integer(n_boot)
  n_partition <- as.integer(n_partition)
  if (!is.logical(return_intermediaries) ||
      length(return_intermediaries) != 1L) {
    stop("'return_intermediaries' must be a single logical value.")
  }
  if (is.null(observables) && n_partition > 1L) {
    message("Coercing 'n_partition' to 1: split-half partitions only vary ",
            "when 'observables' (items mode) are supplied.")
    n_partition <- 1L
  }
  the_sum_fxn <- .lpmec_resolve_partition_aggregation(
    partition_aggregation,
    partition_aggregation_probs
  )
  unit_chr <- if (is.null(unit)) {
    as.character(seq_len(n_obs_input))
  } else {
    as.character(unit)
  }

  subset_measure_input <- function(input, idx) {
    if (is.null(input)) {
      return(NULL)
    }
    if (is.matrix(input) || is.data.frame(input)) {
      return(input[idx, , drop = FALSE])
    }
    if (is.list(input)) {
      return(lapply(input, subset_measure_input, idx = idx))
    }
    input[idx]
  }

  computed <- .lpmec_with_local_seed(seed, {
    runs <- list()
    boot_ids <- integer(0L)
    partition_ids <- integer(0L)
    boot_n_failed <- 0L
    for (boot_i in seq_len(n_boot + 1L)) {
      if (boot_i == 1L) {
        boot_indices <- seq_len(n_obs_input)
        unit_run <- unit_chr
      } else {
        resample <- .lpmec_resample_clusters(unit_chr)
        boot_indices <- resample$indices
        unit_run <- resample$pseudo_unit
      }
      for (partition_i in seq_len(n_partition)) {
        message(sprintf(
          "{boot_i %s of %s} -- {partition_i %s of %s}",
          boot_i, n_boot + 1L, partition_i, n_partition
        ))
        run_onerun <- function() {
          lpmec_panel_onerun(
            Y = if (is.null(Y)) NULL else Y[boot_indices],
            unit = unit_run,
            time = if (is.null(time)) NULL else time[boot_indices],
            observables = subset_measure_input(observables, boot_indices),
            scores = subset_measure_input(scores, boot_indices),
            split_scores = subset_measure_input(split_scores, boot_indices),
            covariates = if (is.null(covariates)) {
              NULL
            } else {
              covariates[boot_indices, , drop = FALSE]
            },
            Y_split_scores = subset_measure_input(Y_split_scores,
                                                  boot_indices),
            design = design,
            diff_k = diff_k,
            estimation_method = estimation_method,
            min_reliability = min_reliability,
            min_cor_n = min_cor_n,
            demean_iterations = demean_iterations,
            ...
          )
        }
        if (boot_i == 1L) {
          # original sample: errors and warnings must surface
          run <- run_onerun()
        } else {
          run <- try(suppressWarnings(run_onerun()), silent = TRUE)
          if (inherits(run, "try-error")) {
            boot_n_failed <- boot_n_failed + 1L
            run <- NULL
          }
        }
        runs[[length(runs) + 1L]] <- run
        boot_ids <- c(boot_ids, boot_i)
        partition_ids <- c(partition_ids, partition_i)
      }
    }
    list(runs = runs, boot_ids = boot_ids, partition_ids = partition_ids,
         boot_n_failed = boot_n_failed)
  })

  runs <- computed$runs
  boot_ids <- computed$boot_ids
  partition_ids <- computed$partition_ids
  template <- runs[[1L]]

  aggregation_fields <- c(
    "ols_coef",
    "corrected_ols_coef", "corrected_ols_coef_split",
    "corrected_ols_coef_triad", "corrected_ols_coef_pair",
    "corrected_ols_lower", "corrected_ols_upper",
    "sd_design_x", "ols_coef_local",
    "corrected_ols_coef_local", "corrected_ols_coef_local_split",
    "corrected_ols_coef_local_triad", "corrected_ols_coef_local_pair",
    "corrected_ols_local_lower", "corrected_ols_local_upper",
    "iv_coef", "corrected_iv_coef", "corrected_iv_coef_within",
    "corrected_iv_coef_cross",
    "split_correlation", "design_split_correlation",
    "split_rho_score", "design_split_rho_score",
    "triad", "design_triad",
    "pair_correlation", "design_pair_correlation",
    "rho_Y_half", "rho_Y",
    "first_stage_fstat",
    "cross_iv_coef", "corrected_cross_iv_coef", "corrected_cross_iv_pair"
  )
  aggregated <- aggregation_fields[
    vapply(aggregation_fields, function(field) {
      length(template[[field]]) > 0L
    }, logical(1L))
  ]

  # NA-fill failed bootstrap runs so per-run matrices stay aligned
  na_template <- lapply(template[aggregated], function(values) {
    values[] <- NA_real_
    values
  })
  failed_runs <- vapply(runs, is.null, logical(1L))
  for (i in which(failed_runs)) {
    runs[[i]] <- na_template
  }
  failed_boots <- unique(boot_ids[failed_runs])

  results <- list(
    measure_names = template$measure_names,
    n_measures = template$n_measures,
    measure_source = template$measure_source,
    has_splits = template$has_splits,
    sign_flipped = template$sign_flipped,
    covariate_names = template$covariate_names,
    design = design,
    diff_k = template$diff_k,
    n_obs = template$n_obs,
    n_units = template$n_units,
    n_boot = n_boot,
    n_partition = n_partition,
    boot_n_failed = computed$boot_n_failed
  )

  intermediaries <- list(
    Intermediary_BootIndex = boot_ids,
    Intermediary_PartitionIndex = partition_ids
  )
  for (field in aggregated) {
    runs_matrix <- .lpmec_runs_matrix(runs, field)
    by_boot <- .lpmec_panel_aggregate_by_boot(runs_matrix, boot_ids,
                                              the_sum_fxn)
    kept_rows <- !(sort(unique(boot_ids)) %in% failed_boots)
    kept_rows[1L] <- TRUE
    by_boot_kept <- by_boot[kept_rows, , drop = FALSE]
    results[[field]] <- stats::setNames(
      as.numeric(by_boot[1L, ]), colnames(by_boot)
    )
    results[[paste0(field, "_se")]] <- .lpmec_boot_sd(by_boot_kept, n_boot)
    results[[paste0(field, "_lower")]] <-
      .lpmec_boot_quantile(by_boot_kept, n_boot, 0.025)
    results[[paste0(field, "_upper")]] <-
      .lpmec_boot_quantile(by_boot_kept, n_boot, 0.975)
    intermediaries[[paste0("Intermediary_", field)]] <- runs_matrix
  }
  for (field in setdiff(aggregation_fields, aggregated)) {
    results[[field]] <- template[[field]]
  }

  results <- c(results, list(
    ols_cluster_se = template$ols_se,
    corrected_ols_source = template$corrected_ols_source,
    corrected_ols_local_source = template$corrected_ols_local_source,
    corrected_iv_source = template$corrected_iv_source,
    reliability = template$reliability,
    cor_matrices = template$cor_matrices,
    cor_n_matrices = template$cor_n_matrices,
    x_est = template$x_est,
    x_est1 = template$x_est1,
    x_est2 = template$x_est2,
    min_reliability = template$min_reliability,
    min_cor_n = template$min_cor_n,
    demean_iterations = template$demean_iterations
  ))
  if (return_intermediaries) {
    results <- c(results, intermediaries)
  }
  class(results) <- "lpmec_panel"
  results
}
