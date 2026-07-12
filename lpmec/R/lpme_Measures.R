# lpme_Measures.R -- measure resolution (items / scores / split scores),
# sign alignment, Spearman-Brown, triad reliabilities, reliability tables,
# Prop-7 sensitivity ranges, and the exported lpmec_reliability_bounds().

#' Pooled z-score with NA-safe moments
#'
#' Non-finite inputs are treated as missing; if the standard deviation is not
#' finite or is zero, an all-\code{NA} vector is returned.
#'
#' @noRd
.lpmec_zscore <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(x * NA_real_)
  }
  (x - mean(x, na.rm = TRUE)) / s
}

#' Resolve the tri-source measure inputs into per-measure score matrices
#'
#' Merges the three measure sources by name into per-measure
#' \code{x_est}/\code{x_est1}/\code{x_est2} columns (pooled z-scored and, by
#' default, sign-aligned):
#' \itemize{
#'   \item \code{observables}: list of item matrices; each measure is scored
#'     by \code{\link{lpmec_onerun}} and its \code{x_est}, \code{x_est1},
#'     \code{x_est2} are harvested (the \code{lpmec_multivariate_onerun}
#'     precedent). When \code{Y} is \code{NULL}, a deterministic placeholder
#'     outcome is used because only the latent scores are harvested.
#'   \item \code{split_scores}: named list of n x 2 half-score matrices; each
#'     half is z-scored, and \code{x_est} is the z-score of the mean of the
#'     two z-scored halves.
#'   \item \code{scores}: named list (or matrix) of full measure scores;
#'     z-scored, with \code{x_est1}/\code{x_est2} left \code{NA}.
#' }
#' A measure name appearing in both \code{scores} and \code{split_scores} is
#' merged (the supplied full score is used as \code{x_est}; the halves supply
#' \code{x_est1}/\code{x_est2}). Names in \code{observables} must not appear
#' in the other sources; duplicates within a source are an error.
#'
#' @return List with \code{measure_names}, \code{n_measures}, \code{n_obs},
#'   \code{x_est}, \code{x_est1}, \code{x_est2} (n x M matrices),
#'   \code{source}, \code{has_splits}, \code{sign_flipped}, and
#'   \code{scalar_runs}.
#'
#' @noRd
.lpmec_resolve_measures <- function(observables = NULL,
                                    scores = NULL,
                                    split_scores = NULL,
                                    Y = NULL,
                                    estimation_method = "averaging",
                                    align_signs = TRUE,
                                    partitions = NULL,
                                    ...) {
  if (is.null(observables) && is.null(scores) && is.null(split_scores)) {
    stop("At least one of 'observables', 'scores', or 'split_scores' is required.")
  }
  if (!is.null(partitions) && !is.list(partitions)) {
    stop("'partitions' must be NULL or a list keyed by measure name.")
  }

  default_names <- function(input_names, n_inputs, prefix) {
    if (is.null(input_names)) {
      input_names <- rep("", n_inputs)
    }
    blank <- !nzchar(input_names) | is.na(input_names)
    input_names[blank] <- paste0(prefix, which(blank))
    input_names
  }

  obs_list <- NULL
  if (!is.null(observables)) {
    obs_list <- .lpmec_as_observables_list(observables)
    names(obs_list) <- default_names(names(obs_list), length(obs_list), "X")
    for (j in seq_along(obs_list)) {
      if (!is.data.frame(obs_list[[j]]) && !is.matrix(obs_list[[j]])) {
        stop("'observables' measure '", names(obs_list)[j],
             "' must be a data.frame or matrix.")
      }
      if (ncol(obs_list[[j]]) < 4L) {
        stop("'observables' measure '", names(obs_list)[j],
             "' has fewer than 4 columns; provide it via 'split_scores' instead.")
      }
    }
  }

  score_list <- NULL
  if (!is.null(scores)) {
    if (is.matrix(scores) || is.data.frame(scores)) {
      score_matrix <- .lpmec_numeric_observable_matrix(scores, label = "'scores'")
      score_list <- lapply(seq_len(ncol(score_matrix)), function(j) score_matrix[, j])
      names(score_list) <- colnames(score_matrix)
    } else if (is.list(scores)) {
      score_list <- scores
    } else if (is.numeric(scores)) {
      score_list <- list(scores)
    } else {
      stop("'scores' must be a numeric vector, matrix, data.frame, or list.")
    }
    names(score_list) <- default_names(names(score_list), length(score_list), "S")
    for (j in seq_along(score_list)) {
      score_j <- score_list[[j]]
      if (is.matrix(score_j) || is.data.frame(score_j)) {
        if (ncol(score_j) != 1L) {
          stop("'scores' measure '", names(score_list)[j],
               "' must be a numeric vector or single-column matrix.")
        }
        score_j <- .lpmec_numeric_observable_matrix(
          score_j,
          label = paste0("'scores' measure '", names(score_list)[j], "'")
        )[, 1L]
      }
      if (!is.numeric(score_j)) {
        stop("'scores' measure '", names(score_list)[j], "' must be numeric.")
      }
      score_j <- as.numeric(score_j)
      score_j[!is.finite(score_j)] <- NA_real_
      score_list[[j]] <- score_j
    }
  }

  split_list <- NULL
  if (!is.null(split_scores)) {
    if (is.matrix(split_scores) || is.data.frame(split_scores)) {
      split_list <- list(split_scores)
    } else if (is.list(split_scores)) {
      split_list <- split_scores
    } else {
      stop("'split_scores' must be a matrix, data.frame, or list of n x 2 matrices.")
    }
    names(split_list) <- default_names(names(split_list), length(split_list), "H")
    for (j in seq_along(split_list)) {
      if (!is.data.frame(split_list[[j]]) && !is.matrix(split_list[[j]])) {
        stop("'split_scores' measure '", names(split_list)[j],
             "' must be a data.frame or matrix.")
      }
      half_matrix <- .lpmec_numeric_observable_matrix(
        split_list[[j]],
        label = paste0("'split_scores' measure '", names(split_list)[j], "'")
      )
      if (ncol(half_matrix) != 2L) {
        stop("'split_scores' measure '", names(split_list)[j],
             "' must have exactly 2 columns (one per half score).")
      }
      split_list[[j]] <- half_matrix
    }
  }

  for (source_name in c("observables", "scores", "split_scores")) {
    source_names <- switch(source_name,
                           observables = names(obs_list),
                           scores = names(score_list),
                           split_scores = names(split_list))
    if (anyDuplicated(source_names) > 0L) {
      stop("Duplicate measure name(s) within '", source_name, "': ",
           paste(unique(source_names[duplicated(source_names)]),
                 collapse = ", "), ".")
    }
  }
  # A measure may appear in BOTH 'scores' and 'split_scores': the full score
  # and its half scores are merged by name (e.g., an official composite index
  # plus component half-scores). 'observables' names must be unique to that
  # source because item scoring already fills all three score columns.
  obs_conflicts <- intersect(names(obs_list),
                             c(names(score_list), names(split_list)))
  if (length(obs_conflicts) > 0L) {
    stop("Measure name(s) supplied in 'observables' cannot also appear in ",
         "'scores' or 'split_scores': ",
         paste(obs_conflicts, collapse = ", "), ".")
  }
  measure_names <- unique(c(names(obs_list), names(score_list),
                            names(split_list)))

  row_counts <- c(
    vapply(obs_list, nrow, integer(1L)),
    vapply(score_list, length, integer(1L)),
    vapply(split_list, nrow, integer(1L))
  )
  if (length(unique(row_counts)) != 1L) {
    stop("All measure inputs must have the same number of rows/observations. ",
         "Received: ", paste(unique(row_counts), collapse = ", "), ".")
  }
  n_obs <- as.integer(row_counts[[1L]])

  n_measures <- length(measure_names)
  x_est <- x_est1 <- x_est2 <- matrix(
    NA_real_, nrow = n_obs, ncol = n_measures,
    dimnames = list(NULL, measure_names)
  )
  scalar_runs <- stats::setNames(vector("list", n_measures), measure_names)
  source <- stats::setNames(rep(NA_character_, n_measures), measure_names)

  if (!is.null(obs_list)) {
    if (is.null(Y)) {
      Y_for_scoring <- as.numeric(seq_len(n_obs))
    } else {
      if (!is.numeric(Y) || length(Y) != n_obs) {
        stop("'Y' must be a numeric vector with one value per row of the ",
             "measure inputs.")
      }
      Y_for_scoring <- as.numeric(Y)
    }
    method_list <- .lpmec_recycle_to_q(
      estimation_method, length(obs_list), "estimation_method"
    )
    for (j in seq_along(obs_list)) {
      measure_name <- names(obs_list)[j]
      run <- lpmec_onerun(
        Y = Y_for_scoring,
        observables = obs_list[[j]],
        estimation_method = method_list[[j]],
        partition = if (is.null(partitions)) NULL else partitions[[measure_name]],
        ...
      )
      x_est[, measure_name] <- .lpmec_zscore(as.numeric(run$x_est))
      x_est1[, measure_name] <- .lpmec_zscore(as.numeric(run$x_est1))
      x_est2[, measure_name] <- .lpmec_zscore(as.numeric(run$x_est2))
      scalar_runs[[measure_name]] <- run
      source[[measure_name]] <- "observables"
    }
  }

  if (!is.null(score_list)) {
    for (j in seq_along(score_list)) {
      measure_name <- names(score_list)[j]
      x_est[, measure_name] <- .lpmec_zscore(score_list[[j]])
      source[[measure_name]] <- "scores"
    }
  }

  if (!is.null(split_list)) {
    for (j in seq_along(split_list)) {
      measure_name <- names(split_list)[j]
      half1 <- .lpmec_zscore(split_list[[j]][, 1L])
      half2 <- .lpmec_zscore(split_list[[j]][, 2L])
      x_est1[, measure_name] <- half1
      x_est2[, measure_name] <- half2
      if (identical(source[[measure_name]], "scores")) {
        # Merged measure: keep the user-supplied full score as x_est and add
        # the half scores for split-based reliabilities.
        source[[measure_name]] <- "scores+split_scores"
      } else {
        x_est[, measure_name] <- .lpmec_zscore(rowMeans(cbind(half1, half2)))
        source[[measure_name]] <- "split_scores"
      }
    }
  }

  if (align_signs) {
    aligned <- .lpmec_align_signs(x_est, x_est1, x_est2)
    x_est <- aligned$x_est
    x_est1 <- aligned$x_est1
    x_est2 <- aligned$x_est2
    sign_flipped <- aligned$flipped
  } else {
    sign_flipped <- stats::setNames(rep(FALSE, n_measures), measure_names)
  }

  list(
    measure_names = measure_names,
    n_measures = n_measures,
    n_obs = n_obs,
    x_est = x_est,
    x_est1 = x_est1,
    x_est2 = x_est2,
    source = source,
    has_splits = source != "scores",
    sign_flipped = sign_flipped,
    scalar_runs = scalar_runs
  )
}

#' Sign-align measure score columns before any correlations are taken
#'
#' Each measure's \code{x_est} column is flipped when it correlates negatively
#' (pairwise-complete) with the first informative measure column, and the two
#' half-score columns are then flipped to correlate positively with their own
#' measure's \code{x_est}. Columns whose reference correlation is not finite
#' are left unchanged.
#'
#' @return List with aligned \code{x_est}, \code{x_est1}, \code{x_est2}, and
#'   a named logical vector \code{flipped}.
#'
#' @noRd
.lpmec_align_signs <- function(x_est, x_est1 = NULL, x_est2 = NULL) {
  x_est <- as.matrix(x_est)
  if (!is.null(x_est1)) {
    x_est1 <- as.matrix(x_est1)
  }
  if (!is.null(x_est2)) {
    x_est2 <- as.matrix(x_est2)
  }
  n_measures <- ncol(x_est)
  flipped <- stats::setNames(rep(FALSE, n_measures), colnames(x_est))

  pairwise_cor <- function(a, b) {
    suppressWarnings(stats::cor(a, b, use = "pairwise.complete.obs"))
  }

  informative <- which(vapply(seq_len(n_measures), function(j) {
    sum(is.finite(x_est[, j])) >= 2L
  }, logical(1L)))
  if (length(informative) >= 1L && n_measures >= 2L) {
    reference <- x_est[, informative[1L]]
    for (j in setdiff(seq_len(n_measures), informative[1L])) {
      r <- pairwise_cor(x_est[, j], reference)
      if (is.finite(r) && r < 0) {
        x_est[, j] <- -x_est[, j]
        if (!is.null(x_est1)) {
          x_est1[, j] <- -x_est1[, j]
        }
        if (!is.null(x_est2)) {
          x_est2[, j] <- -x_est2[, j]
        }
        flipped[j] <- TRUE
      }
    }
  }

  for (j in seq_len(n_measures)) {
    if (!is.null(x_est1)) {
      r1 <- pairwise_cor(x_est1[, j], x_est[, j])
      if (is.finite(r1) && r1 < 0) {
        x_est1[, j] <- -x_est1[, j]
      }
    }
    if (!is.null(x_est2)) {
      r2 <- pairwise_cor(x_est2[, j], x_est[, j])
      if (is.finite(r2) && r2 < 0) {
        x_est2[, j] <- -x_est2[, j]
      }
    }
  }

  list(x_est = x_est, x_est1 = x_est1, x_est2 = x_est2, flipped = flipped)
}

#' Spearman-Brown step-up of a (mean) split correlation
#'
#' \code{rho = factor * r / (1 + (factor - 1) * r)}; with \code{factor = 2}
#' this steps a half-score reliability up to the full-score reliability.
#' Non-finite inputs propagate to \code{NA}.
#'
#' @noRd
.lpmec_spearman_brown <- function(r, factor = 2) {
  if (!is.numeric(factor) || length(factor) != 1L || !is.finite(factor) ||
      factor < 1) {
    stop("'factor' must be a single numeric value greater than or equal to 1.")
  }
  r <- as.numeric(r)
  out <- (factor * r) / (1 + (factor - 1) * r)
  out[!is.finite(r)] <- NA_real_
  out
}

#' Triad ("triangulation") reliabilities from a cross-measure correlation matrix
#'
#' For measure m and two other measures l, k (Prop 7b):
#' \code{rho*_m = r_ml * r_mk / r_lk}, which is consistent when the measures'
#' total errors are pairwise orthogonal. Requires at least 3 measures;
#' otherwise all entries are \code{NA}. With more than 3 measures the finite
#' candidates over all pairs \{l, k\} are averaged.
#'
#' @noRd
.lpmec_triad_reliabilities <- function(cor_matrix) {
  cor_matrix <- as.matrix(cor_matrix)
  n_measures <- ncol(cor_matrix)
  measure_names <- colnames(cor_matrix)
  if (is.null(measure_names)) {
    measure_names <- paste0("X", seq_len(n_measures))
  }
  out <- stats::setNames(rep(NA_real_, n_measures), measure_names)
  if (n_measures < 3L) {
    return(out)
  }
  for (m in seq_len(n_measures)) {
    others <- setdiff(seq_len(n_measures), m)
    other_pairs <- utils::combn(others, 2L)
    candidates <- apply(other_pairs, 2L, function(pair) {
      cor_matrix[m, pair[1L]] * cor_matrix[m, pair[2L]] /
        cor_matrix[pair[1L], pair[2L]]
    })
    candidates <- candidates[is.finite(candidates)]
    if (length(candidates) > 0L) {
      out[m] <- mean(candidates)
    }
  }
  out
}

#' Long reliability table over designs and measures
#'
#' For each design transform: transforms all measure score columns, computes
#' gated within-measure split correlations, their Spearman-Brown step-up to
#' the full-score scale (\code{rho_split}), the gated cross-measure
#' correlation matrix, and the per-measure triad reliabilities. The raw
#' split correlation estimates the reliability of a \emph{half} score;
#' \code{rho_split = 2r/(1+r)} places it on the scale of the full score that
#' the corrections divide by (exact when the full score is the average of
#' two independent parallel halves; a parallel-forms approximation
#' otherwise).
#'
#' @param measures Output of \code{.lpmec_resolve_measures}.
#'
#' @return List with \code{table} (data.frame: design, measure, source,
#'   split_correlation, split_n, rho_split, triad), \code{cor_matrices}, and
#'   \code{cor_n_matrices} (named by design).
#'
#' @noRd
.lpmec_reliability_table <- function(measures,
                                     unit,
                                     time = NULL,
                                     designs = "pooled",
                                     diff_k = 1L,
                                     demean_iterations = 25L,
                                     min_cor_n = 30L) {
  measure_names <- measures$measure_names
  n_measures <- measures$n_measures
  rows <- list()
  cor_matrices <- list()
  cor_n_matrices <- list()

  for (design_name in designs) {
    transformed <- .lpmec_panel_transform(
      measures$x_est, unit = unit, time = time, design = design_name,
      diff_k = diff_k, demean_iterations = demean_iterations
    )
    transformed1 <- .lpmec_panel_transform(
      measures$x_est1, unit = unit, time = time, design = design_name,
      diff_k = diff_k, demean_iterations = demean_iterations
    )
    transformed2 <- .lpmec_panel_transform(
      measures$x_est2, unit = unit, time = time, design = design_name,
      diff_k = diff_k, demean_iterations = demean_iterations
    )

    split_correlation <- rep(NA_real_, n_measures)
    split_n <- rep(NA_integer_, n_measures)
    for (m in seq_len(n_measures)) {
      split_stats <- .lpmec_cor_transformed(
        transformed1[, m], transformed2[, m], unit = unit,
        design = "pooled", min_cor_n = min_cor_n
      )
      split_correlation[m] <- split_stats[["r"]]
      split_n[m] <- as.integer(split_stats[["n"]])
    }

    cross_matrix <- diag(1, n_measures)
    dimnames(cross_matrix) <- list(measure_names, measure_names)
    cross_n <- matrix(NA_integer_, n_measures, n_measures,
                      dimnames = list(measure_names, measure_names))
    diag(cross_n) <- colSums(is.finite(transformed))
    if (n_measures >= 2L) {
      for (l in seq_len(n_measures - 1L)) {
        for (m in seq.int(l + 1L, n_measures)) {
          cross_stats <- .lpmec_cor_transformed(
            transformed[, l], transformed[, m], unit = unit,
            design = "pooled", min_cor_n = min_cor_n
          )
          cross_matrix[l, m] <- cross_matrix[m, l] <- cross_stats[["r"]]
          cross_n[l, m] <- cross_n[m, l] <- as.integer(cross_stats[["n"]])
        }
      }
    }
    triads <- .lpmec_triad_reliabilities(cross_matrix)

    rows[[design_name]] <- data.frame(
      design = design_name,
      measure = measure_names,
      source = unname(measures$source),
      split_correlation = split_correlation,
      split_n = split_n,
      rho_split = as.numeric(.lpmec_spearman_brown(split_correlation)),
      triad = unname(triads),
      stringsAsFactors = FALSE
    )
    cor_matrices[[design_name]] <- cross_matrix
    cor_n_matrices[[design_name]] <- cross_n
  }

  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  list(
    table = table,
    cor_matrices = cor_matrices,
    cor_n_matrices = cor_n_matrices
  )
}

#' Prop-7 sensitivity range from a reliability table
#'
#' Per (measure, design) row, \code{rho_lo}/\code{rho_hi} are the minimum
#' and maximum of the finite reliability candidates \{triad,
#' score-scale split reliability \code{rho_split}\}. Both candidates
#' estimate the full score's reliability, so the range compares
#' like-for-like; it is a sensitivity range, not an identified set, unless
#' the corresponding orthogonality or ratio condition of Proposition 7 is
#' maintained. Candidates below \code{min_reliability} (including negative
#' values) are excluded -- those corrections are not identified in
#' practice -- and a single warning lists the affected rows.
#'
#' @return The input table with \code{rho_lo} and \code{rho_hi} appended.
#'
#' @noRd
.lpmec_bounds_from_reliability <- function(reliability_table,
                                           min_reliability = 0.05,
                                           warn = TRUE) {
  n_rows <- nrow(reliability_table)
  rho_lo <- rho_hi <- rep(NA_real_, n_rows)
  floored_labels <- character(0L)
  rho_split_column <- if (is.null(reliability_table$rho_split)) {
    as.numeric(.lpmec_spearman_brown(reliability_table$split_correlation))
  } else {
    reliability_table$rho_split
  }
  for (i in seq_len(n_rows)) {
    candidates <- c(reliability_table$triad[i],
                    rho_split_column[i])
    finite_candidates <- candidates[is.finite(candidates)]
    usable <- finite_candidates[finite_candidates >= min_reliability]
    if (length(usable) < length(finite_candidates)) {
      floored_labels <- c(
        floored_labels,
        paste0(reliability_table$measure[i],
               " (", reliability_table$design[i], ")")
      )
    }
    if (length(usable) > 0L) {
      rho_lo[i] <- min(usable)
      rho_hi[i] <- max(usable)
    }
  }
  if (warn && length(floored_labels) > 0L) {
    warning(
      "Reliability estimate(s) below 'min_reliability' (", min_reliability,
      ") were set to NA when forming bounds for: ",
      paste(unique(floored_labels), collapse = ", "), ".",
      call. = FALSE
    )
  }
  reliability_table$rho_lo <- rho_lo
  reliability_table$rho_hi <- rho_hi
  reliability_table
}

#' Reliability sensitivity diagnostics for latent measures under panel
#' design transforms
#'
#' Estimates, for each latent measure and each requested design transform,
#' the within-measure split correlation (with its Spearman-Brown step-up to
#' the full-score scale, \code{rho_split}) and the cross-measure triad
#' reliability, and reports the sensitivity range
#' \code{[rho_lo, rho_hi]} spanned by the two candidates for the measure's
#' design reliability (Proposition 7 of the accompanying working paper).
#' The range is a sensitivity diagnostic, not a partial-identification
#' interval, unless Proposition 7's corresponding orthogonality or ratio
#' condition is maintained. Optionally attaches a cluster (unit) bootstrap
#' for all reported reliability quantities.
#'
#' @param observables Optional list of item matrices or data frames, one per
#'   measure (a single matrix is treated as one measure). Each measure is
#'   scored via \code{\link{lpmec_onerun}} with \code{estimation_method}, and
#'   its full-battery and split-half latent scores are harvested. Each item
#'   matrix must have at least 4 columns; measures with fewer components must
#'   be supplied through \code{split_scores}.
#' @param scores Optional named list (or matrix/data frame with one column
#'   per measure) of pre-computed full measure scores. Scores are pooled
#'   z-scored; split-based quantities are unavailable for these measures.
#' @param split_scores Optional named list of n x 2 matrices of half scores,
#'   one per measure. Each half is pooled z-scored and the full-measure score
#'   is the z-score of the mean of the two z-scored halves. This is the
#'   required input mode for measures with only 2 components.
#' @param unit Optional vector of unit (cluster) identifiers, one per row.
#'   Required for the \code{"within"}, \code{"twoway"}, and
#'   \code{"difference"} designs. When \code{NULL}, each row is its own
#'   cluster and only \code{"pooled"} designs are allowed.
#' @param time Optional numeric vector of time identifiers, one per row.
#'   Required for the \code{"twoway"} and \code{"difference"} designs.
#'   (unit, time) pairs must be unique.
#' @param designs Character vector of design transforms; any of
#'   \code{"pooled"} (identity), \code{"within"} (unit demeaning),
#'   \code{"twoway"} (iterated unit-and-time demeaning), and
#'   \code{"difference"} (exact-gap k-period differencing). Default
#'   \code{"pooled"}.
#' @param diff_k Positive integer gap for the \code{"difference"} design.
#'   Rows lacking a same-unit observation at exactly \code{time - diff_k} are
#'   set to \code{NA} rather than differenced against a nearer observation.
#' @param estimation_method Estimation method(s) passed to
#'   \code{\link{lpmec_onerun}} for \code{observables}-mode measures (single
#'   value or one per measure). Default \code{"averaging"}.
#' @param min_reliability Reliability floor. Estimated reliabilities below
#'   this value (or negative) are excluded from the bounds and reported as
#'   \code{NA} with a single warning, since corrections that divide by such
#'   values are not identified in practice. Default \code{0.05}.
#' @param min_cor_n Minimum number of complete pairs required to report any
#'   correlation; correlations on fewer pairs are \code{NA}. Default
#'   \code{30}.
#' @param demean_iterations Number of iterated demeaning passes for the
#'   \code{"twoway"} design (one pass is exact only on balanced panels).
#'   Default \code{25}.
#' @param n_boot Non-negative integer number of cluster-bootstrap
#'   replications. Units are resampled with replacement and receive fresh
#'   pseudo-ids (a unit drawn twice enters as two distinct clusters); the
#'   whole pipeline (scoring, sign alignment, transforms, reliabilities,
#'   bounds) is re-run on each draw. Default \code{0}.
#' @param seed Optional seed applied locally to the scoring partitions and
#'   the bootstrap (the caller's random-number state is restored on exit).
#' @param ... Additional arguments passed to \code{\link{lpmec_onerun}} for
#'   \code{observables}-mode measures (e.g., \code{ordinal},
#'   \code{mcmc_control}).
#'
#' @return A list of class \code{lpmec_reliability_bounds} containing:
#' \itemize{
#'   \item \code{reliability}: data frame with one row per (design, measure)
#'     holding \code{split_correlation}, \code{split_n}, \code{rho_split}
#'     (the Spearman-Brown step-up of the split correlation to the
#'     full-score scale), \code{triad}, \code{rho_lo}, \code{rho_hi}, and --
#'     when \code{n_boot > 0} -- bootstrap \code{_se}, \code{_lower}, and
#'     \code{_upper} columns for each of these reliability quantities.
#'   \item \code{cor_matrices} and \code{cor_n_matrices}: per-design
#'     cross-measure correlation matrices and complete-pair counts.
#'   \item \code{x_est}, \code{x_est1}, \code{x_est2}: pooled z-scored,
#'     sign-aligned measure score matrices (halves are \code{NA} for
#'     \code{scores}-mode measures).
#'   \item \code{measure_names}, \code{n_measures}, \code{designs},
#'     \code{measure_source}, \code{sign_flipped}, \code{n_obs},
#'     \code{n_units}, \code{diff_k}, \code{demean_iterations},
#'     \code{min_reliability}, \code{min_cor_n}, \code{n_boot},
#'     \code{n_boot_failed}.
#'   \item When \code{n_boot > 0}: \code{Intermediary_BootIndex} and
#'     per-replication matrices \code{Intermediary_split_correlation},
#'     \code{Intermediary_rho_split}, \code{Intermediary_triad},
#'     \code{Intermediary_rho_lo}, \code{Intermediary_rho_hi} (row 1 is the
#'     original sample).
#' }
#'
#' @details
#' The split correlation between two half scores of the same measure
#' consistently estimates the reliability of a \emph{half} score;
#' \code{rho_split = 2r/(1+r)} steps it up to the scale of the full score
#' (exact when the full score is the average of two independent parallel
#' halves, a parallel-forms approximation otherwise). The two candidates
#' entering the range therefore both target the full score's reliability.
#' Their directional interpretation is model-dependent (Proposition 7):
#' when the halves share a common systematic error component orthogonal to
#' the trait, the split-based value exceeds construct-relevant reliability
#' (Prop 7a); the triad \code{rho*_m = r_ml * r_mk / r_lk} equals it exactly
#' when the measures' total errors are pairwise orthogonal (Prop 7b) and is
#' a lower bound only under the target-specific ratio condition of Prop 7c
#' -- shared error among the other two measures pushes it down, but shared
#' error involving the target measure can push it up. Absent those
#' conditions, \code{[rho_lo, rho_hi]} is the [min, max] of two sensitivity
#' candidates, not an identified set. With fewer than 3 measures the triad
#' is not available and both endpoints collapse to \code{rho_split}.
#' Correlations are computed after pooled z-scoring and sign alignment of
#' all measure scores; correlations with fewer than \code{min_cor_n}
#' complete pairs, and reliabilities below \code{min_reliability}, are
#' reported as \code{NA} rather than propagated into unstable divisions.
#'
#' @examples
#' \donttest{
#' set.seed(100)
#' n <- 400
#' latent <- rnorm(n)
#' half1 <- latent + rnorm(n, sd = 0.6)
#' half2 <- latent + rnorm(n, sd = 0.6)
#' score2 <- latent + rnorm(n, sd = 0.7)
#' score3 <- latent + rnorm(n, sd = 0.7)
#'
#' bounds <- lpmec_reliability_bounds(
#'   split_scores = list(measure1 = cbind(half1, half2)),
#'   scores = list(measure2 = score2, measure3 = score3)
#' )
#' bounds$reliability
#' }
#'
#' @export
lpmec_reliability_bounds <- function(observables = NULL,
                                     scores = NULL,
                                     split_scores = NULL,
                                     unit = NULL,
                                     time = NULL,
                                     designs = "pooled",
                                     diff_k = 1L,
                                     estimation_method = "averaging",
                                     min_reliability = 0.05,
                                     min_cor_n = 30L,
                                     demean_iterations = 25L,
                                     n_boot = 0L,
                                     seed = NULL,
                                     ...) {
  if (is.null(observables) && is.null(scores) && is.null(split_scores)) {
    stop("At least one of 'observables', 'scores', or 'split_scores' is required.")
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || !is.finite(n_boot) ||
      n_boot != floor(n_boot) || n_boot < 0) {
    stop("'n_boot' must be a single non-negative integer.")
  }
  n_boot <- as.integer(n_boot)

  input_rows <- function(input) {
    if (is.null(input)) {
      return(NULL)
    }
    if (is.matrix(input) || is.data.frame(input)) {
      return(nrow(input))
    }
    if (is.list(input)) {
      if (length(input) < 1L) {
        return(NULL)
      }
      return(input_rows(input[[1L]]))
    }
    length(input)
  }
  n_obs <- c(input_rows(observables), input_rows(scores),
             input_rows(split_scores))[1L]
  if (is.null(n_obs) || !is.finite(n_obs) || n_obs < 1L) {
    stop("Could not determine the number of observations from the measure inputs.")
  }

  prep <- .lpmec_prepare_panel_inputs(
    n_obs = n_obs,
    unit = unit,
    time = time,
    design = designs,
    diff_k = diff_k,
    min_reliability = min_reliability,
    min_cor_n = min_cor_n,
    demean_iterations = demean_iterations
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

  run_pipeline <- function(observables_run, scores_run, split_scores_run,
                           unit_run, time_run) {
    measures <- .lpmec_resolve_measures(
      observables = observables_run,
      scores = scores_run,
      split_scores = split_scores_run,
      Y = NULL,
      estimation_method = estimation_method,
      ...
    )
    reliability <- .lpmec_reliability_table(
      measures = measures,
      unit = unit_run,
      time = time_run,
      designs = prep$design,
      diff_k = prep$diff_k,
      demean_iterations = prep$demean_iterations,
      min_cor_n = prep$min_cor_n
    )
    list(measures = measures, reliability = reliability)
  }

  computed <- .lpmec_with_local_seed(seed, {
    point <- run_pipeline(observables, scores, split_scores,
                          prep$unit, prep$time)
    point_table <- .lpmec_bounds_from_reliability(
      point$reliability$table,
      min_reliability = prep$min_reliability,
      warn = TRUE
    )
    boot_tables <- vector("list", n_boot)
    for (boot_i in seq_len(n_boot)) {
      resample <- .lpmec_resample_clusters(prep$unit)
      boot_table <- try(suppressWarnings({
        boot_run <- run_pipeline(
          subset_measure_input(observables, resample$indices),
          subset_measure_input(scores, resample$indices),
          subset_measure_input(split_scores, resample$indices),
          resample$pseudo_unit,
          prep$time[resample$indices]
        )
        .lpmec_bounds_from_reliability(
          boot_run$reliability$table,
          min_reliability = prep$min_reliability,
          warn = FALSE
        )
      }), silent = TRUE)
      if (!inherits(boot_table, "try-error")) {
        boot_tables[[boot_i]] <- boot_table
      }
    }
    list(point = point, point_table = point_table, boot_tables = boot_tables)
  })

  point <- computed$point
  point_table <- computed$point_table
  boot_tables <- computed$boot_tables

  stat_fields <- c("split_correlation", "rho_split", "triad",
                   "rho_lo", "rho_hi")
  intermediaries <- NULL
  if (n_boot > 0L) {
    column_keys <- paste0(point_table$measure, ".", point_table$design)
    intermediaries <- lapply(stat_fields, function(field) {
      values <- matrix(NA_real_, nrow = n_boot + 1L,
                       ncol = nrow(point_table),
                       dimnames = list(NULL, column_keys))
      values[1L, ] <- point_table[[field]]
      for (boot_i in seq_len(n_boot)) {
        if (!is.null(boot_tables[[boot_i]])) {
          values[boot_i + 1L, ] <- boot_tables[[boot_i]][[field]]
        }
      }
      values
    })
    names(intermediaries) <- stat_fields
    for (field in stat_fields) {
      values <- intermediaries[[field]]
      point_table[[paste0(field, "_se")]] <-
        as.numeric(.lpmec_boot_sd(values, n_boot))
      point_table[[paste0(field, "_lower")]] <-
        as.numeric(.lpmec_boot_quantile(values, n_boot, 0.025))
      point_table[[paste0(field, "_upper")]] <-
        as.numeric(.lpmec_boot_quantile(values, n_boot, 0.975))
    }
  }

  results <- list(
    measure_names = point$measures$measure_names,
    n_measures = point$measures$n_measures,
    designs = prep$design,
    reliability = point_table,
    cor_matrices = point$reliability$cor_matrices,
    cor_n_matrices = point$reliability$cor_n_matrices,
    x_est = point$measures$x_est,
    x_est1 = point$measures$x_est1,
    x_est2 = point$measures$x_est2,
    measure_source = point$measures$source,
    sign_flipped = point$measures$sign_flipped,
    n_obs = prep$n_obs,
    n_units = length(unique(prep$unit)),
    diff_k = prep$diff_k,
    demean_iterations = prep$demean_iterations,
    min_reliability = prep$min_reliability,
    min_cor_n = prep$min_cor_n,
    n_boot = n_boot,
    n_boot_failed = if (n_boot > 0L) {
      sum(vapply(boot_tables, is.null, logical(1L)))
    } else {
      0L
    }
  )
  if (n_boot > 0L) {
    results$Intermediary_BootIndex <- seq_len(n_boot + 1L)
    results$Intermediary_split_correlation <- intermediaries$split_correlation
    results$Intermediary_rho_split <- intermediaries$rho_split
    results$Intermediary_triad <- intermediaries$triad
    results$Intermediary_rho_lo <- intermediaries$rho_lo
    results$Intermediary_rho_hi <- intermediaries$rho_hi
  }
  class(results) <- "lpmec_reliability_bounds"
  results
}
