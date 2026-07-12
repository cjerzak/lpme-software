# lpme_PanelTransforms.R -- internal panel/design-transform helpers shared by
# lpmec_reliability_bounds() and the lpmec_panel_*() family. Pure base R +
# stats; cluster-robust covariances via sandwich::vcovCL and split-IV fits via
# AER::ivreg. Ported from the verified V2 reference implementations
# (apply_transform / cor_transformed / resample_countries).

.lpmec_valid_panel_designs <- function() {
  c("pooled", "within", "twoway", "difference")
}

#' Validate and normalize shared panel inputs
#'
#' Checks the unit/time identifiers, requested design transform(s), and the
#' shared tuning constants used by the panel and reliability-bounds functions.
#' When \code{unit} is \code{NULL}, only pooled designs are allowed and each
#' row is treated as its own cluster (so the cluster bootstrap reduces to a
#' row bootstrap).
#'
#' @param n_obs Single positive integer number of rows.
#' @param unit Optional vector of unit (cluster) identifiers of length
#'   \code{n_obs}.
#' @param time Optional numeric vector of time identifiers of length
#'   \code{n_obs}; required for the \code{"twoway"} and \code{"difference"}
#'   designs.
#' @param design Character vector of design transforms; each element must be
#'   one of \code{"pooled"}, \code{"within"}, \code{"twoway"},
#'   \code{"difference"}.
#' @param diff_k Single positive integer gap for the difference design.
#' @param min_reliability Single numeric floor in [0, 1) below which
#'   reliability estimates are treated as unidentified.
#' @param min_cor_n Single positive integer minimum number of complete pairs
#'   required to report a correlation.
#' @param demean_iterations Single positive integer number of iterated
#'   two-way demeaning passes.
#'
#' @return A list with normalized \code{n_obs}, \code{unit} (character),
#'   \code{time} (numeric or \code{NULL}), \code{design}, \code{diff_k},
#'   \code{min_reliability}, \code{min_cor_n}, \code{demean_iterations}, and
#'   \code{unit_supplied}.
#'
#' @noRd
.lpmec_prepare_panel_inputs <- function(n_obs,
                                        unit = NULL,
                                        time = NULL,
                                        design = "pooled",
                                        diff_k = 1L,
                                        min_reliability = 0.05,
                                        min_cor_n = 30L,
                                        demean_iterations = 25L) {
  if (!is.numeric(n_obs) || length(n_obs) != 1L || !is.finite(n_obs) ||
      n_obs != floor(n_obs) || n_obs < 1L) {
    stop("'n_obs' must be a single positive integer.")
  }
  n_obs <- as.integer(n_obs)

  valid_designs <- .lpmec_valid_panel_designs()
  if (!is.character(design) || length(design) < 1L || anyNA(design)) {
    stop("'design' must be a character vector of design transforms.")
  }
  design <- unique(design)
  unknown_designs <- setdiff(design, valid_designs)
  if (length(unknown_designs) > 0L) {
    stop("'design' must contain only: ", paste(valid_designs, collapse = ", "),
         ". Received: '", paste(unknown_designs, collapse = "', '"), "'")
  }

  unit_supplied <- !is.null(unit)
  if (unit_supplied) {
    if (length(unit) != n_obs) {
      stop("'unit' must have length ", n_obs, ". Received: ", length(unit), ".")
    }
    if (anyNA(unit)) {
      stop("'unit' must not contain missing values.")
    }
    unit <- as.character(unit)
  } else {
    non_pooled <- setdiff(design, "pooled")
    if (length(non_pooled) > 0L) {
      stop("'unit' is required for design(s): ",
           paste(non_pooled, collapse = ", "), ".")
    }
    unit <- as.character(seq_len(n_obs))
  }

  if (!is.null(time)) {
    if (!is.numeric(time)) {
      stop("'time' must be a numeric vector.")
    }
    if (length(time) != n_obs) {
      stop("'time' must have length ", n_obs, ". Received: ", length(time), ".")
    }
    if (anyNA(time)) {
      stop("'time' must not contain missing values.")
    }
    time <- as.numeric(time)
    if (anyDuplicated(paste(unit, time, sep = "\r")) > 0L) {
      stop("Duplicate (unit, time) pairs found; each unit-time combination ",
           "must appear at most once.")
    }
  } else {
    needs_time <- intersect(design, c("twoway", "difference"))
    if (length(needs_time) > 0L) {
      stop("'time' is required for design(s): ",
           paste(needs_time, collapse = ", "), ".")
    }
  }

  if (!is.numeric(diff_k) || length(diff_k) != 1L || !is.finite(diff_k) ||
      diff_k != floor(diff_k) || diff_k < 1L) {
    stop("'diff_k' must be a single positive integer.")
  }
  if (!is.numeric(min_reliability) || length(min_reliability) != 1L ||
      !is.finite(min_reliability) || min_reliability < 0 ||
      min_reliability >= 1) {
    stop("'min_reliability' must be a single finite value in [0, 1).")
  }
  if (!is.numeric(min_cor_n) || length(min_cor_n) != 1L ||
      !is.finite(min_cor_n) || min_cor_n != floor(min_cor_n) ||
      min_cor_n < 1L) {
    stop("'min_cor_n' must be a single positive integer.")
  }
  if (!is.numeric(demean_iterations) || length(demean_iterations) != 1L ||
      !is.finite(demean_iterations) ||
      demean_iterations != floor(demean_iterations) ||
      demean_iterations < 1L) {
    stop("'demean_iterations' must be a single positive integer.")
  }

  list(
    n_obs = n_obs,
    unit = unit,
    time = time,
    design = design,
    diff_k = as.integer(diff_k),
    min_reliability = as.numeric(min_reliability),
    min_cor_n = as.integer(min_cor_n),
    demean_iterations = as.integer(demean_iterations),
    unit_supplied = unit_supplied
  )
}

#' NA-safe group demeaning of each column of a numeric matrix
#'
#' Subtracts group means computed with \code{na.rm = TRUE}, so rows with
#' missing values do not contaminate the group means and remain missing.
#'
#' @noRd
.lpmec_group_demean <- function(x_mat, groups) {
  for (j in seq_len(ncol(x_mat))) {
    group_means <- stats::ave(x_mat[, j], groups,
                              FUN = function(v) mean(v, na.rm = TRUE))
    x_mat[, j] <- x_mat[, j] - group_means
  }
  x_mat
}

#' Apply a panel design transform to a vector or matrix of columns
#'
#' Designs: \code{"pooled"} is the identity (location shifts are irrelevant
#' for the downstream correlations); \code{"within"} demeans by unit with
#' NA-safe group means; \code{"twoway"} applies iterated unit-then-time
#' demeaning (\code{demean_iterations} passes; one pass is exact only on
#' balanced panels); \code{"difference"} takes the exact-gap k-difference:
#' each row is matched to the row of the same unit at exactly
#' \code{time - diff_k}, and rows lacking an exact partner become \code{NA}.
#' Non-finite results are normalized to \code{NA}.
#'
#' @param x Numeric vector or matrix (rows are observations).
#' @param unit Unit identifiers (length \code{nrow(x)}).
#' @param time Numeric time identifiers; required for \code{"twoway"} and
#'   \code{"difference"}.
#' @param design One design transform.
#' @param diff_k Positive integer gap for \code{"difference"}.
#' @param demean_iterations Number of iterated two-way demeaning passes.
#'
#' @return Object of the same shape as \code{x} (vector in, vector out).
#'
#' @noRd
.lpmec_panel_transform <- function(x,
                                   unit,
                                   time = NULL,
                                   design = c("pooled", "within", "twoway", "difference"),
                                   diff_k = 1L,
                                   demean_iterations = 25L) {
  design <- match.arg(design)
  is_vector_input <- is.null(dim(x))
  x_mat <- as.matrix(x)
  storage.mode(x_mat) <- "double"
  if (nrow(x_mat) != length(unit)) {
    stop("'x' must have one row per element of 'unit'.")
  }
  unit <- as.character(unit)
  if (design %in% c("twoway", "difference")) {
    if (is.null(time)) {
      stop("'time' is required for design = '", design, "'.")
    }
    if (length(time) != nrow(x_mat)) {
      stop("'time' must have one value per row of 'x'.")
    }
    time <- as.numeric(time)
  }

  if (design == "within") {
    x_mat <- .lpmec_group_demean(x_mat, unit)
  } else if (design == "twoway") {
    for (iteration in seq_len(demean_iterations)) {
      x_mat <- .lpmec_group_demean(x_mat, unit)
      x_mat <- .lpmec_group_demean(x_mat, time)
    }
  } else if (design == "difference") {
    row_key <- paste(unit, time, sep = "\r")
    partner_key <- paste(unit, time - diff_k, sep = "\r")
    partner_row <- match(partner_key, row_key)
    x_mat <- x_mat - x_mat[partner_row, , drop = FALSE]
  }

  x_mat[!is.finite(x_mat)] <- NA_real_
  if (is_vector_input) {
    return(as.numeric(x_mat))
  }
  x_mat
}

#' Correlation of two columns under a design transform, with a sample gate
#'
#' Both columns are transformed jointly, then correlated over complete pairs.
#' When fewer than \code{min_cor_n} complete pairs remain, the correlation is
#' reported as \code{NA} (the pair count is always reported).
#'
#' @return Named numeric vector \code{c(r = ..., n = ...)}.
#'
#' @noRd
.lpmec_cor_transformed <- function(a,
                                   b,
                                   unit,
                                   time = NULL,
                                   design = "pooled",
                                   diff_k = 1L,
                                   demean_iterations = 25L,
                                   min_cor_n = 30L) {
  transformed <- .lpmec_panel_transform(
    cbind(a = as.numeric(a), b = as.numeric(b)),
    unit = unit,
    time = time,
    design = design,
    diff_k = diff_k,
    demean_iterations = demean_iterations
  )
  complete <- is.finite(transformed[, 1L]) & is.finite(transformed[, 2L])
  n_pairs <- sum(complete)
  if (n_pairs < min_cor_n) {
    return(c(r = NA_real_, n = n_pairs))
  }
  r <- suppressWarnings(
    stats::cor(transformed[complete, 1L], transformed[complete, 2L])
  )
  if (!is.finite(r)) {
    r <- NA_real_
  }
  c(r = r, n = n_pairs)
}

#' Cluster-robust OLS on (already transformed) panel data
#'
#' Fits \code{stats::lm} of \code{Y} on \code{x} (plus optional covariates)
#' over complete cases and reports the latent slope with a cluster-robust
#' standard error from \code{sandwich::vcovCL}. \code{x_sd} is the sample
#' standard deviation of \code{x} on the estimation sample, used to place
#' coefficients on the design-local scale (Proposition 2a).
#'
#' @return List with \code{coef}, \code{se}, \code{coef_all}, \code{se_all},
#'   \code{n_obs}, \code{n_clusters}, \code{x_sd}.
#'
#' @noRd
.lpmec_panel_ols <- function(Y, x, covariates = NULL, cluster) {
  Y <- as.numeric(Y)
  x <- as.numeric(x)
  n <- length(Y)
  if (length(x) != n) {
    stop("'x' must have the same length as 'Y'.")
  }
  if (length(cluster) != n) {
    stop("'cluster' must have the same length as 'Y'.")
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

  complete <- is.finite(Y) & is.finite(x)
  if (ncol(covariate_matrix) > 0L) {
    complete <- complete & rowSums(!is.finite(covariate_matrix)) == 0L
  }
  n_complete <- sum(complete)
  x_sd <- if (n_complete >= 2L) stats::sd(x[complete]) else NA_real_
  empty_result <- list(
    coef = NA_real_,
    se = NA_real_,
    coef_all = NULL,
    se_all = NULL,
    n_obs = n_complete,
    n_clusters = length(unique(cluster[complete])),
    x_sd = x_sd
  )
  if (n_complete < ncol(covariate_matrix) + 3L) {
    return(empty_result)
  }

  fit_data <- as.data.frame(
    cbind(Y = Y[complete], x = x[complete],
          covariate_matrix[complete, , drop = FALSE])
  )
  colnames(fit_data) <- make.names(colnames(fit_data), unique = TRUE)
  fit <- try(stats::lm(Y ~ ., data = fit_data), silent = TRUE)
  if (inherits(fit, "try-error") || anyNA(stats::coef(fit))) {
    return(empty_result)
  }
  vcov_cluster <- try(
    sandwich::vcovCL(fit, cluster = factor(cluster[complete])),
    silent = TRUE
  )
  if (inherits(vcov_cluster, "try-error")) {
    return(empty_result)
  }
  coef_all <- stats::coef(fit)
  se_all <- sqrt(diag(vcov_cluster))
  list(
    coef = unname(coef_all[["x"]]),
    se = unname(se_all[["x"]]),
    coef_all = coef_all[-1L],
    se_all = se_all[-1L],
    n_obs = n_complete,
    n_clusters = length(unique(cluster[complete])),
    x_sd = x_sd
  )
}

#' Cluster-robust just-identified IV on (already transformed) panel data
#'
#' Fits \code{AER::ivreg} of \code{Y} on \code{x} instrumented by \code{z}
#' (optional covariates enter both stages), with cluster-robust standard
#' errors from \code{sandwich::vcovCL}. The clustered first-stage F statistic
#' is computed as \code{(coef / se_cl)^2} from the \code{stats::lm} first
#' stage with \code{sandwich::vcovCL}.
#'
#' @return List with \code{coef}, \code{se}, \code{first_stage_fstat},
#'   \code{n_obs}, \code{n_clusters}.
#'
#' @noRd
.lpmec_panel_iv <- function(Y, x, z, covariates = NULL, cluster) {
  Y <- as.numeric(Y)
  x <- as.numeric(x)
  z <- as.numeric(z)
  n <- length(Y)
  if (length(x) != n || length(z) != n) {
    stop("'x' and 'z' must have the same length as 'Y'.")
  }
  if (length(cluster) != n) {
    stop("'cluster' must have the same length as 'Y'.")
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

  complete <- is.finite(Y) & is.finite(x) & is.finite(z)
  if (ncol(covariate_matrix) > 0L) {
    complete <- complete & rowSums(!is.finite(covariate_matrix)) == 0L
  }
  n_complete <- sum(complete)
  empty_result <- list(
    coef = NA_real_,
    se = NA_real_,
    first_stage_fstat = NA_real_,
    n_obs = n_complete,
    n_clusters = length(unique(cluster[complete]))
  )
  if (n_complete < ncol(covariate_matrix) + 4L) {
    return(empty_result)
  }

  fit_data <- as.data.frame(
    cbind(Y = Y[complete], x = x[complete], z = z[complete],
          covariate_matrix[complete, , drop = FALSE])
  )
  colnames(fit_data) <- make.names(colnames(fit_data), unique = TRUE)
  covariate_names <- setdiff(colnames(fit_data), c("Y", "x", "z"))
  rhs <- paste(c("x", covariate_names), collapse = " + ")
  instruments <- paste(c("z", covariate_names), collapse = " + ")
  iv_formula <- stats::as.formula(paste("Y ~", rhs, "|", instruments))
  cluster_factor <- factor(cluster[complete])

  fit <- try(AER::ivreg(iv_formula, data = fit_data), silent = TRUE)
  if (inherits(fit, "try-error") || anyNA(stats::coef(fit))) {
    return(empty_result)
  }
  vcov_cluster <- try(
    sandwich::vcovCL(fit, cluster = cluster_factor),
    silent = TRUE
  )
  if (inherits(vcov_cluster, "try-error")) {
    return(empty_result)
  }

  first_stage_formula <- stats::as.formula(
    paste("x ~", paste(c("z", covariate_names), collapse = " + "))
  )
  first_stage_fstat <- NA_real_
  first_stage <- try(stats::lm(first_stage_formula, data = fit_data),
                     silent = TRUE)
  if (!inherits(first_stage, "try-error") &&
      !anyNA(stats::coef(first_stage))) {
    first_stage_vcov <- try(
      sandwich::vcovCL(first_stage, cluster = cluster_factor),
      silent = TRUE
    )
    if (!inherits(first_stage_vcov, "try-error")) {
      first_stage_coef <- stats::coef(first_stage)[["z"]]
      first_stage_se <- sqrt(diag(first_stage_vcov))[["z"]]
      if (is.finite(first_stage_se) && first_stage_se > 0) {
        first_stage_fstat <- (first_stage_coef / first_stage_se)^2
      }
    }
  }

  list(
    coef = unname(stats::coef(fit)[["x"]]),
    se = unname(sqrt(diag(vcov_cluster))[["x"]]),
    first_stage_fstat = first_stage_fstat,
    n_obs = n_complete,
    n_clusters = length(unique(cluster[complete]))
  )
}

#' Cluster (unit) bootstrap resample with fresh pseudo-ids
#'
#' Samples units with replacement and assigns each draw a fresh pseudo-id, so
#' a unit drawn twice enters as two distinct clusters/fixed-effect groups
#' (port of the V2 \code{resample_countries} convention).
#'
#' @param unit Vector of unit identifiers (one per row).
#'
#' @return List with \code{indices} (row indices of the resampled data, in
#'   draw order) and \code{pseudo_unit} (fresh pseudo-id per resampled row).
#'
#' @noRd
.lpmec_resample_clusters <- function(unit) {
  if (length(unit) < 1L) {
    stop("'unit' must be a non-empty vector.")
  }
  if (anyNA(unit)) {
    stop("'unit' must not contain missing values.")
  }
  unit <- as.character(unit)
  rows_by_unit <- split(seq_along(unit), unit)
  drawn_units <- sample(names(rows_by_unit), length(rows_by_unit),
                        replace = TRUE)
  drawn_rows <- lapply(drawn_units, function(u) rows_by_unit[[u]])
  list(
    indices = unlist(drawn_rows, use.names = FALSE),
    pseudo_unit = rep(
      paste0(drawn_units, "_", seq_along(drawn_units)),
      lengths(drawn_rows)
    )
  )
}
