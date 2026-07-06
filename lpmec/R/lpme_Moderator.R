# lpme_Moderator.R -- latent-moderator interaction correction for randomized
# treatments: exported lpmec_moderator_onerun() and lpmec_moderator(), plus the
# internal normal-equations interaction fit .lpmec_fit_interaction() (port of
# the verified V2 fit_interaction reference implementation).

#' Validate the treatment vector for the moderator functions
#'
#' The treatment must be a numeric vector with one finite value per element of
#' \code{Y} and positive variance. Treatments with more than 2 unique values
#' are allowed (they are treated as continuous) but trigger a message when
#' \code{notify = TRUE}.
#'
#' @noRd
.lpmec_validate_treatment <- function(treatment, n_obs, notify = TRUE) {
  if (is.null(treatment)) {
    stop("'treatment' is required and cannot be NULL.")
  }
  if (!is.numeric(treatment)) {
    stop("'treatment' must be a numeric vector.")
  }
  if (length(treatment) != n_obs) {
    stop("'treatment' must have the same length as 'Y'. Received: ",
         length(treatment), " versus ", n_obs, ".")
  }
  if (any(!is.finite(treatment))) {
    stop("'treatment' must contain only finite values.")
  }
  treatment <- as.numeric(treatment)
  if (stats::var(treatment) == 0) {
    stop("'treatment' must have positive variance (a constant treatment ",
         "cannot identify treatment-effect moderation).")
  }
  if (notify && length(unique(treatment)) > 2L) {
    message("'treatment' has more than 2 unique values; ",
            "treating it as a continuous treatment.")
  }
  invisible(treatment)
}

#' Normal-equations OLS of Y ~ 1 + treatment + x + treatment:x (+ covariates)
#'
#' Port of the V2 \code{fit_interaction} reference: the interaction design is
#' built explicitly and fit with the package's normal-equations OLS
#' (\code{.lpmec_fit_ols}), whose coefficients and classical standard errors
#' match \code{stats::lm(Y ~ treatment * x + ...)}. Rows with non-finite
#' \code{Y}, \code{treatment}, \code{x}, or covariate values are dropped.
#'
#' @param Y Numeric outcome vector.
#' @param treatment Numeric treatment vector.
#' @param x Numeric moderator-score vector.
#' @param covariates Optional numeric covariate matrix (already prepared; no
#'   intercept column).
#'
#' @return List with \code{treatment_coef}/\code{treatment_se},
#'   \code{main_coef}/\code{main_se} (the moderator main effect),
#'   \code{interaction_coef}/\code{interaction_se}, \code{coef_all},
#'   \code{se_all} (all named, intercept omitted), and \code{n_obs} (complete
#'   cases used).
#'
#' @noRd
.lpmec_fit_interaction <- function(Y, treatment, x, covariates = NULL) {
  Y <- as.numeric(Y)
  treatment <- as.numeric(treatment)
  x <- as.numeric(x)
  n <- length(Y)
  if (length(treatment) != n) {
    stop("'treatment' must have the same length as 'Y'.")
  }
  if (length(x) != n) {
    stop("'x' must have the same length as 'Y'.")
  }
  covariate_matrix <- if (is.null(covariates)) {
    matrix(nrow = n, ncol = 0L)
  } else {
    as.matrix(covariates)
  }
  storage.mode(covariate_matrix) <- "double"
  if (nrow(covariate_matrix) != n) {
    stop("'covariates' must have the same number of rows as 'Y'.")
  }
  if (ncol(covariate_matrix) > 0L && is.null(colnames(covariate_matrix))) {
    colnames(covariate_matrix) <- paste0("C", seq_len(ncol(covariate_matrix)))
  }

  complete <- is.finite(Y) & is.finite(treatment) & is.finite(x)
  if (ncol(covariate_matrix) > 0L) {
    complete <- complete & rowSums(!is.finite(covariate_matrix)) == 0L
  }
  n_complete <- sum(complete)
  if (n_complete < ncol(covariate_matrix) + 5L) {
    stop("Too few complete observations (", n_complete,
         ") to fit the interaction regression.")
  }

  regressors <- cbind(
    treatment = treatment[complete],
    x = x[complete],
    "treatment:x" = treatment[complete] * x[complete]
  )
  fit <- .lpmec_fit_ols(
    Y[complete],
    regressors,
    covariate_matrix[complete, , drop = FALSE]
  )
  list(
    treatment_coef = unname(fit$latent_coef[1L]),
    treatment_se = unname(fit$latent_se[1L]),
    main_coef = unname(fit$latent_coef[2L]),
    main_se = unname(fit$latent_se[2L]),
    interaction_coef = unname(fit$latent_coef[3L]),
    interaction_se = unname(fit$latent_se[3L]),
    coef_all = fit$coef_all,
    se_all = fit$se_all,
    n_obs = n_complete
  )
}

#' Single-run latent-moderator interaction correction
#'
#' Estimates the interaction between a treatment and a latent moderator that
#' is measured with error, and corrects the naive interaction coefficient for
#' the attenuation induced by the moderator score's reliability. With a
#' randomized treatment and a standardized moderator score \code{x} of
#' reliability \code{rho}, the naive interaction coefficient from
#' \code{Y ~ treatment * x} converges to \code{b_TX * sqrt(rho)}; the
#' corrected estimate divides by \code{sqrt(rho_score)}, where
#' \code{rho_score} is the Spearman-Brown as-used reliability of the score
#' actually entered in the regression.
#'
#' @param Y Numeric outcome vector with only finite values.
#' @param treatment Numeric treatment vector of the same length as \code{Y}
#'   with only finite values and positive variance. Treatments with more than
#'   2 unique values are allowed and treated as continuous (a message is
#'   emitted).
#' @param observables Optional list of item matrices or data frames, one per
#'   measure (a single matrix is treated as one measure). Each measure is
#'   scored via \code{\link{lpmec_onerun}} with \code{estimation_method}, and
#'   its full-battery and split-half latent scores are harvested. Each item
#'   matrix must have at least 4 columns; measures with fewer components must
#'   be supplied through \code{split_scores}.
#' @param scores Optional named list (or matrix/data frame with one column
#'   per measure) of pre-computed full measure scores. Scores are pooled
#'   z-scored; split-based reliabilities are unavailable for these measures.
#' @param split_scores Optional named list of n x 2 matrices of half scores,
#'   one per measure. Each half is pooled z-scored and the full-measure score
#'   is the z-score of the mean of the two z-scored halves.
#' @param covariates Optional matrix or data frame of observed covariates
#'   included additively in the interaction regression.
#' @param estimation_method Estimation method(s) passed to
#'   \code{\link{lpmec_onerun}} for \code{observables}-mode measures (single
#'   value or one per measure). Default \code{"averaging"}.
#' @param min_reliability Reliability floor in [0, 1). When
#'   \code{rho_score} (or a per-measure triad reliability) is below this
#'   value or not finite, the corresponding corrected coefficient is reported
#'   as \code{NA} with a single warning rather than dividing by a tiny or
#'   undefined reliability. Default \code{0.05}.
#' @param ... Additional arguments passed to \code{\link{lpmec_onerun}} for
#'   \code{observables}-mode measures (e.g., \code{ordinal},
#'   \code{mcmc_control}).
#'
#' @return A list of class \code{lpmec_moderator_onerun} containing:
#' \itemize{
#'   \item \code{treatment_coef}/\code{treatment_se},
#'     \code{main_coef}/\code{main_se},
#'     \code{interaction_coef}/\code{interaction_se}: naive OLS coefficients
#'     and classical standard errors from
#'     \code{Y ~ treatment + x + treatment:x} (plus covariates), where
#'     \code{x} is the moderator score in \code{x_used}.
#'   \item \code{coef_all}, \code{se_all}: all named coefficients (intercept
#'     omitted), including covariates.
#'   \item \code{rho_half}: with one measure, the split-half correlation of
#'     that measure's two half scores; with \code{M >= 2} measures, the mean
#'     pairwise correlation among the z-scored, sign-aligned measure scores.
#'   \item \code{sb_factor}, \code{rho_score}: the Spearman-Brown step-up
#'     factor (2 for half scores of one battery, \code{M} for \code{M}
#'     averaged measures) and the as-used score reliability
#'     \code{rho_score = sb_factor * rho_half / (1 + (sb_factor - 1) *
#'     rho_half)}.
#'   \item \code{correction_factor}: \code{sqrt(rho_score)} when usable,
#'     otherwise \code{NA}.
#'   \item \code{corrected_interaction_coef}, \code{corrected_main_coef}:
#'     \code{interaction_coef / sqrt(rho_score)} and
#'     \code{main_coef / sqrt(rho_score)} (\code{NA} when \code{rho_score} is
#'     below \code{min_reliability} or not finite).
#'   \item \code{reliability_floored}: logical, \code{TRUE} when the
#'     headline correction was floored to \code{NA}.
#'   \item \code{per_measure}: when \code{M >= 3}, a data frame with one row
#'     per measure holding the naive per-measure interaction coefficient (and
#'     classical standard error), the pooled triad reliability
#'     \code{rho*_m}, and the triad-corrected interaction coefficient;
#'     \code{NULL} otherwise.
#'   \item \code{cor_matrix}: pooled cross-measure correlation matrix
#'     (\code{NULL} when \code{M == 1}).
#'   \item \code{x_used}: the moderator score entered in the regression (the
#'     z-scored full-battery score for one measure; the z-scored mean of the
#'     aligned z-scored measure scores for \code{M >= 2}).
#'   \item \code{x_est}, \code{x_est1}, \code{x_est2}: per-measure pooled
#'     z-scored, sign-aligned score matrices (halves are \code{NA} for
#'     \code{scores}-mode measures).
#'   \item \code{measure_names}, \code{n_measures}, \code{measure_source},
#'     \code{sign_flipped}, \code{covariate_names}, \code{n_obs},
#'     \code{n_obs_used}, \code{min_reliability}, \code{scalar_runs}.
#' }
#'
#' @details
#' The correction divides by the square root of the reliability of the score
#' actually used as the regressor (Proposition 5 of the accompanying working
#' paper). With a single item battery, the split-half correlation
#' \code{rho_half} estimates the reliability of a half score, so the
#' full-battery score reliability is the Spearman-Brown step-up
#' \code{2 * rho_half / (1 + rho_half)}; dividing by \code{sqrt(rho_half)}
#' instead would overcorrect. With \code{M >= 2} measures, the regressor is
#' the mean of the \code{M} aligned z-scored measure scores and the analogous
#' generalized step-up \code{M * rbar / (1 + (M - 1) * rbar)} is applied to
#' the mean pairwise correlation \code{rbar}. When \code{M >= 3}, the
#' per-measure triad reliabilities \code{rho*_m = r_ml * r_mk / r_lk} --
#' consistent when measurement errors are independent across measures --
#' provide an alternative correction reported in \code{per_measure}.
#' Covariates enter the regression additively; the correction is exact when
#' the covariates are independent of the latent moderator and approximate
#' otherwise. Reliabilities below \code{min_reliability} (or not estimable,
#' as with a single \code{scores}-mode measure) yield \code{NA} corrected
#' coefficients and a single warning.
#'
#' @examples
#' \donttest{
#' set.seed(100)
#' n <- 600
#' X <- rnorm(n)
#' treatment <- rbinom(n, 1, 0.5)
#' Y <- 0.2 * treatment + 0.2 * X + 0.3 * treatment * X + rnorm(n)
#' items <- matrix(rbinom(n * 4L, 1,
#'                        stats::pnorm(outer(X, c(1.2, 1, 0.8, 1.1)) - 0.2)),
#'                 n, 4L)
#' colnames(items) <- paste0("item", 1:4)
#'
#' run <- lpmec_moderator_onerun(
#'   Y = Y,
#'   treatment = treatment,
#'   observables = items,
#'   estimation_method = "averaging"
#' )
#' c(naive = run$interaction_coef,
#'   rho_score = run$rho_score,
#'   corrected = run$corrected_interaction_coef)
#' }
#'
#' @export
lpmec_moderator_onerun <- function(Y,
                                   treatment,
                                   observables = NULL,
                                   scores = NULL,
                                   split_scores = NULL,
                                   covariates = NULL,
                                   estimation_method = "averaging",
                                   min_reliability = 0.05,
                                   ...) {
  if (missing(Y) || is.null(Y)) {
    stop("'Y' is required and cannot be NULL.")
  }
  if (!is.numeric(Y)) {
    stop("'Y' must be a numeric vector.")
  }
  if (length(Y) < 10L) {
    stop("'Y' must have at least 10 observations. Received: ", length(Y))
  }
  if (any(!is.finite(Y))) {
    stop("'Y' must contain only finite values.")
  }
  Y <- as.numeric(Y)
  if (missing(treatment)) {
    stop("'treatment' is required and cannot be NULL.")
  }
  .lpmec_validate_treatment(treatment, length(Y), notify = TRUE)
  treatment <- as.numeric(treatment)
  if (!is.numeric(min_reliability) || length(min_reliability) != 1L ||
      !is.finite(min_reliability) || min_reliability < 0 ||
      min_reliability >= 1) {
    stop("'min_reliability' must be a single finite value in [0, 1).")
  }
  if (is.null(observables) && is.null(scores) && is.null(split_scores)) {
    stop("At least one of 'observables', 'scores', or 'split_scores' is required.")
  }

  covariate_matrix <- .lpmec_prepare_covariates(covariates, length(Y))
  measures <- .lpmec_resolve_measures(
    observables = observables,
    scores = scores,
    split_scores = split_scores,
    Y = Y,
    estimation_method = estimation_method,
    ...
  )
  if (measures$n_obs != length(Y)) {
    stop("All measure inputs must have one row per element of 'Y'. ",
         "Received ", measures$n_obs, " rows for ", length(Y),
         " observations.")
  }
  n_measures <- measures$n_measures

  pairwise_cor <- function(a, b) {
    r <- suppressWarnings(stats::cor(a, b, use = "pairwise.complete.obs"))
    if (!is.finite(r)) NA_real_ else r
  }

  if (n_measures == 1L) {
    x_used <- measures$x_est[, 1L]
    rho_half <- pairwise_cor(measures$x_est1[, 1L], measures$x_est2[, 1L])
    sb_factor <- 2
    cor_matrix <- NULL
  } else {
    x_used <- .lpmec_zscore(rowMeans(measures$x_est))
    cor_matrix <- suppressWarnings(
      stats::cor(measures$x_est, use = "pairwise.complete.obs")
    )
    off_diagonal <- cor_matrix[upper.tri(cor_matrix)]
    off_diagonal <- off_diagonal[is.finite(off_diagonal)]
    rho_half <- if (length(off_diagonal) > 0L) {
      mean(off_diagonal)
    } else {
      NA_real_
    }
    sb_factor <- n_measures
  }
  rho_score <- as.numeric(.lpmec_spearman_brown(rho_half, factor = sb_factor))

  fit <- .lpmec_fit_interaction(Y, treatment, x_used, covariate_matrix)

  floored_labels <- character(0L)
  rho_usable <- is.finite(rho_score) && rho_score >= min_reliability
  if (!rho_usable) {
    floored_labels <- "rho_score"
  }
  correction_factor <- if (rho_usable) sqrt(rho_score) else NA_real_
  corrected_interaction_coef <- if (rho_usable) {
    fit$interaction_coef / correction_factor
  } else {
    NA_real_
  }
  corrected_main_coef <- if (rho_usable) {
    fit$main_coef / correction_factor
  } else {
    NA_real_
  }

  per_measure <- NULL
  if (n_measures >= 3L) {
    triads <- .lpmec_triad_reliabilities(cor_matrix)
    measure_interaction_coef <- rep(NA_real_, n_measures)
    measure_interaction_se <- rep(NA_real_, n_measures)
    measure_corrected <- rep(NA_real_, n_measures)
    for (m in seq_len(n_measures)) {
      fit_m <- try(
        .lpmec_fit_interaction(Y, treatment, measures$x_est[, m],
                               covariate_matrix),
        silent = TRUE
      )
      if (!inherits(fit_m, "try-error")) {
        measure_interaction_coef[m] <- fit_m$interaction_coef
        measure_interaction_se[m] <- fit_m$interaction_se
      }
      triad_usable <- is.finite(triads[m]) && triads[m] >= min_reliability
      if (triad_usable) {
        measure_corrected[m] <- measure_interaction_coef[m] / sqrt(triads[m])
      } else {
        floored_labels <- c(
          floored_labels,
          paste0("triad (", measures$measure_names[m], ")")
        )
      }
    }
    per_measure <- data.frame(
      measure = measures$measure_names,
      interaction_coef = measure_interaction_coef,
      interaction_se = measure_interaction_se,
      triad = unname(triads),
      corrected_interaction_coef = measure_corrected,
      stringsAsFactors = FALSE
    )
  }

  if (length(floored_labels) > 0L) {
    warning(
      "Reliability estimate(s) below 'min_reliability' (", min_reliability,
      ") or not estimable; the corresponding corrected coefficient(s) are ",
      "NA: ", paste(floored_labels, collapse = ", "), ".",
      call. = FALSE
    )
  }

  results <- list(
    n_obs = length(Y),
    n_obs_used = fit$n_obs,
    n_measures = n_measures,
    measure_names = measures$measure_names,
    measure_source = measures$source,
    sign_flipped = measures$sign_flipped,
    covariate_names = colnames(covariate_matrix),
    treatment_coef = fit$treatment_coef,
    treatment_se = fit$treatment_se,
    main_coef = fit$main_coef,
    main_se = fit$main_se,
    interaction_coef = fit$interaction_coef,
    interaction_se = fit$interaction_se,
    coef_all = fit$coef_all,
    se_all = fit$se_all,
    rho_half = rho_half,
    sb_factor = sb_factor,
    rho_score = rho_score,
    correction_factor = correction_factor,
    corrected_interaction_coef = corrected_interaction_coef,
    corrected_main_coef = corrected_main_coef,
    reliability_floored = !rho_usable,
    min_reliability = min_reliability,
    per_measure = per_measure,
    cor_matrix = cor_matrix,
    x_used = x_used,
    x_est = measures$x_est,
    x_est1 = measures$x_est1,
    x_est2 = measures$x_est2,
    scalar_runs = measures$scalar_runs
  )
  class(results) <- "lpmec_moderator_onerun"
  results
}

#' Aggregated latent-moderator interaction correction with bootstrap
#'
#' Runs \code{\link{lpmec_moderator_onerun}} over repeated split-half
#' partitions (in \code{observables} mode) and row (or stratified) bootstrap
#' samples, re-running the whole pipeline -- scoring with fresh splits,
#' reliability estimation, interaction regression, and correction -- on every
#' (bootstrap, partition) draw so the reported uncertainty reflects all
#' estimation steps.
#'
#' @inheritParams lpmec_moderator_onerun
#' @param n_boot Non-negative integer number of n-out-of-n bootstrap
#'   replications. Default \code{32}.
#' @param n_partition Positive integer number of fresh split-half partitions
#'   per original or bootstrap sample. Fresh partitions require item-level
#'   inputs, so when \code{observables} is \code{NULL} the value is coerced
#'   to 1 with a message. Default \code{10}.
#' @param partition_aggregation Aggregation strategy across partitions within
#'   each bootstrap sample; see \code{\link{lpmec}}. Default \code{"median"}.
#' @param partition_aggregation_probs Quantile probabilities for winsorized
#'   or trimmed partition aggregation.
#' @param boot_basis Optional vector of length \code{length(Y)} of strata
#'   labels for the bootstrap. With all-unique values (the default,
#'   \code{seq_along(Y)}), rows are resampled with replacement; otherwise
#'   rows are resampled with replacement within each stratum (for example,
#'   \code{boot_basis = treatment} preserves the arm sizes).
#' @param return_intermediaries Logical. If \code{TRUE}, returns per-run
#'   vectors/matrices of all aggregated quantities.
#' @param seed Optional seed applied locally to the partitions and the
#'   bootstrap (the caller's random-number state is restored on exit).
#'
#' @return A list of class \code{lpmec_moderator} containing, for each of
#'   \code{treatment_coef}, \code{main_coef}, \code{interaction_coef},
#'   \code{corrected_main_coef}, \code{corrected_interaction_coef},
#'   \code{rho_half}, and \code{rho_score}: the original-sample estimate
#'   (aggregated across partitions) plus bootstrap \code{_se},
#'   \code{_lower}, and \code{_upper} summaries (2.5\% and 97.5\%
#'   percentiles) when \code{n_boot >= 1}. Also included are
#'   \code{coef_all}/\code{coef_se_all}, \code{sb_factor},
#'   \code{per_measure} and \code{cor_matrix} from the original-sample run,
#'   \code{x_used}/\code{x_est}/\code{x_est1}/\code{x_est2},
#'   \code{measure_names}, \code{n_measures}, \code{measure_source},
#'   \code{covariate_names}, \code{n_obs}, \code{n_boot},
#'   \code{n_partition}, \code{min_reliability}, and -- when
#'   \code{return_intermediaries = TRUE} -- \code{Intermediary_BootIndex},
#'   \code{Intermediary_PartitionIndex}, and per-run
#'   \code{Intermediary_*} vectors (plus \code{Intermediary_coef_all}).
#'
#' @details
#' The first run (bootstrap index 1) is the original sample; its
#' partition-aggregated estimates are the reported point estimates, and the
#' remaining \code{n_boot} runs feed the bootstrap standard errors and
#' percentile intervals. Because reliability estimation is re-done on every
#' draw, the intervals for \code{corrected_interaction_coef} propagate the
#' sampling uncertainty in \code{rho_half} and \code{rho_score} -- both of
#' which are reported with their own bootstrap summaries as headline
#' quantities. Warnings and messages from the pipeline are shown only for
#' the first (original-sample, first-partition) run; corrected coefficients
#' floored to \code{NA} in bootstrap replications propagate \code{NA}
#' through the default median aggregation.
#'
#' @examples
#' \donttest{
#' set.seed(101)
#' n <- 400
#' X <- rnorm(n)
#' treatment <- rbinom(n, 1, 0.5)
#' Y <- 0.2 * treatment + 0.2 * X + 0.3 * treatment * X + rnorm(n)
#' items <- matrix(rbinom(n * 4L, 1,
#'                        stats::pnorm(outer(X, c(1.2, 1, 0.8, 1.1)) - 0.2)),
#'                 n, 4L)
#' colnames(items) <- paste0("item", 1:4)
#'
#' result <- lpmec_moderator(
#'   Y = Y,
#'   treatment = treatment,
#'   observables = items,
#'   estimation_method = "averaging",
#'   n_boot = 4L,
#'   n_partition = 2L,
#'   seed = 7
#' )
#' c(corrected = result$corrected_interaction_coef,
#'   lower = result$corrected_interaction_lower,
#'   upper = result$corrected_interaction_upper)
#' }
#'
#' @export
lpmec_moderator <- function(Y,
                            treatment,
                            observables = NULL,
                            scores = NULL,
                            split_scores = NULL,
                            covariates = NULL,
                            n_boot = 32L,
                            n_partition = 10L,
                            partition_aggregation = "median",
                            partition_aggregation_probs = c(0.01, 0.99),
                            boot_basis = seq_along(Y),
                            return_intermediaries = TRUE,
                            estimation_method = "averaging",
                            min_reliability = 0.05,
                            seed = NULL,
                            ...) {
  if (missing(Y) || is.null(Y)) {
    stop("'Y' is required and cannot be NULL.")
  }
  if (!is.numeric(Y)) {
    stop("'Y' must be a numeric vector.")
  }
  if (missing(treatment)) {
    stop("'treatment' is required and cannot be NULL.")
  }
  .lpmec_validate_treatment(treatment, length(Y), notify = FALSE)
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
  if (length(boot_basis) != length(Y)) {
    stop("'boot_basis' must have the same length as 'Y'.")
  }
  if (!is.logical(return_intermediaries) || length(return_intermediaries) != 1L) {
    stop("'return_intermediaries' must be a single logical value.")
  }
  if (!is.null(covariates) && !is.matrix(covariates) &&
      !is.data.frame(covariates)) {
    stop("'covariates' must be NULL, a data.frame, or a matrix.")
  }

  if (is.null(observables) && n_partition > 1L) {
    message("'n_partition' applies only when 'observables' are supplied ",
            "(fresh split-half partitions require item-level inputs); ",
            "coercing 'n_partition' to 1.")
    n_partition <- 1L
  }

  the_sum_fxn <- .lpmec_resolve_partition_aggregation(
    partition_aggregation,
    partition_aggregation_probs
  )

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

  scalar_fields <- c(
    "treatment_coef", "main_coef", "interaction_coef",
    "corrected_main_coef", "corrected_interaction_coef",
    "rho_half", "rho_score"
  )

  computed <- .lpmec_with_local_seed(seed, {
    runs <- list()
    boot_ids <- integer(0L)
    partition_ids <- integer(0L)

    for (boot_i in seq_len(n_boot + 1L)) {
      if (boot_i == 1L) {
        boot_indices <- seq_along(Y)
      } else if (length(unique(boot_basis)) == length(Y)) {
        boot_indices <- sample(seq_along(Y), length(Y), replace = TRUE)
      } else {
        strata <- split(seq_along(boot_basis), as.character(boot_basis))
        boot_indices <- unlist(
          lapply(strata, function(idx) sample(idx, length(idx), replace = TRUE)),
          use.names = FALSE
        )
      }

      for (partition_i in seq_len(n_partition)) {
        message(sprintf(
          "{boot_i %s of %s} -- {partition_i %s of %s}",
          boot_i,
          n_boot + 1L,
          partition_i,
          n_partition
        ))

        run_once <- function() {
          lpmec_moderator_onerun(
            Y = Y[boot_indices],
            treatment = treatment[boot_indices],
            observables = subset_measure_input(observables, boot_indices),
            scores = subset_measure_input(scores, boot_indices),
            split_scores = subset_measure_input(split_scores, boot_indices),
            covariates = if (is.null(covariates)) {
              NULL
            } else {
              covariates[boot_indices, , drop = FALSE]
            },
            estimation_method = estimation_method,
            min_reliability = min_reliability,
            ...
          )
        }
        run <- if (boot_i == 1L && partition_i == 1L) {
          run_once()
        } else {
          suppressWarnings(suppressMessages(run_once()))
        }
        run$scalar_summary <- vapply(
          scalar_fields,
          function(field) as.numeric(run[[field]]),
          numeric(1L)
        )
        runs[[length(runs) + 1L]] <- run
        boot_ids <- c(boot_ids, boot_i)
        partition_ids <- c(partition_ids, partition_i)
      }
    }
    list(runs = runs, boot_ids = boot_ids, partition_ids = partition_ids)
  })

  runs <- computed$runs
  boot_ids <- computed$boot_ids
  partition_ids <- computed$partition_ids

  scalar_runs_matrix <- .lpmec_runs_matrix(runs, "scalar_summary")
  coef_all_runs <- .lpmec_runs_matrix(runs, "coef_all")

  scalar_by_boot <- .lpmec_aggregate_by_boot(
    scalar_runs_matrix, boot_ids, the_sum_fxn
  )
  coef_all_by_boot <- .lpmec_aggregate_by_boot(
    coef_all_runs, boot_ids, the_sum_fxn
  )

  point <- scalar_by_boot[1L, ]
  boot_se <- .lpmec_boot_sd(scalar_by_boot, n_boot)
  boot_lower <- .lpmec_boot_quantile(scalar_by_boot, n_boot, 0.025)
  boot_upper <- .lpmec_boot_quantile(scalar_by_boot, n_boot, 0.975)

  results <- list(
    n_obs = runs[[1L]]$n_obs,
    n_measures = runs[[1L]]$n_measures,
    measure_names = runs[[1L]]$measure_names,
    measure_source = runs[[1L]]$measure_source,
    sign_flipped = runs[[1L]]$sign_flipped,
    covariate_names = runs[[1L]]$covariate_names,
    treatment_coef = unname(point[["treatment_coef"]]),
    treatment_se = unname(boot_se[["treatment_coef"]]),
    treatment_lower = unname(boot_lower[["treatment_coef"]]),
    treatment_upper = unname(boot_upper[["treatment_coef"]]),
    main_coef = unname(point[["main_coef"]]),
    main_se = unname(boot_se[["main_coef"]]),
    main_lower = unname(boot_lower[["main_coef"]]),
    main_upper = unname(boot_upper[["main_coef"]]),
    interaction_coef = unname(point[["interaction_coef"]]),
    interaction_se = unname(boot_se[["interaction_coef"]]),
    interaction_lower = unname(boot_lower[["interaction_coef"]]),
    interaction_upper = unname(boot_upper[["interaction_coef"]]),
    corrected_main_coef = unname(point[["corrected_main_coef"]]),
    corrected_main_se = unname(boot_se[["corrected_main_coef"]]),
    corrected_main_lower = unname(boot_lower[["corrected_main_coef"]]),
    corrected_main_upper = unname(boot_upper[["corrected_main_coef"]]),
    corrected_interaction_coef = unname(point[["corrected_interaction_coef"]]),
    corrected_interaction_se = unname(boot_se[["corrected_interaction_coef"]]),
    corrected_interaction_lower = unname(boot_lower[["corrected_interaction_coef"]]),
    corrected_interaction_upper = unname(boot_upper[["corrected_interaction_coef"]]),
    rho_half = unname(point[["rho_half"]]),
    rho_half_se = unname(boot_se[["rho_half"]]),
    rho_half_lower = unname(boot_lower[["rho_half"]]),
    rho_half_upper = unname(boot_upper[["rho_half"]]),
    rho_score = unname(point[["rho_score"]]),
    rho_score_se = unname(boot_se[["rho_score"]]),
    rho_score_lower = unname(boot_lower[["rho_score"]]),
    rho_score_upper = unname(boot_upper[["rho_score"]]),
    sb_factor = runs[[1L]]$sb_factor,
    coef_all = coef_all_by_boot[1L, ],
    coef_se_all = .lpmec_boot_sd(coef_all_by_boot, n_boot),
    coef_lower_all = .lpmec_boot_quantile(coef_all_by_boot, n_boot, 0.025),
    coef_upper_all = .lpmec_boot_quantile(coef_all_by_boot, n_boot, 0.975),
    per_measure = runs[[1L]]$per_measure,
    cor_matrix = runs[[1L]]$cor_matrix,
    x_used = runs[[1L]]$x_used,
    x_est = runs[[1L]]$x_est,
    x_est1 = runs[[1L]]$x_est1,
    x_est2 = runs[[1L]]$x_est2,
    n_boot = n_boot,
    n_partition = n_partition,
    min_reliability = min_reliability
  )

  if (return_intermediaries) {
    intermediaries <- list(
      Intermediary_BootIndex = boot_ids,
      Intermediary_PartitionIndex = partition_ids
    )
    for (field in scalar_fields) {
      intermediaries[[paste0("Intermediary_", field)]] <-
        unname(scalar_runs_matrix[, field])
    }
    intermediaries$Intermediary_coef_all <- coef_all_runs
    results <- c(results, intermediaries)
  }

  class(results) <- "lpmec_moderator"
  results
}
