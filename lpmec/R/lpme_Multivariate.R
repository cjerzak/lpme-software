#' Multivariate latent-predictor measurement-error correction
#'
#' Implements the split-indicator multivariate IV correction for several latent
#' predictors. Each latent predictor is estimated from its own indicator matrix.
#'
#' @param Y Numeric outcome vector.
#' @param observables A list of matrices or data frames, one per latent predictor.
#' @param covariates Optional matrix or data frame of observed covariates.
#' @param observables_groupings Optional list of grouping vectors, one per latent
#'   predictor. Defaults to the column names of each observable matrix.
#' @param make_observables_groupings Logical scalar or vector passed to
#'   \code{\link{lpmec_onerun}} for each latent predictor.
#' @param estimation_method Character scalar or vector passed to
#'   \code{\link{lpmec_onerun}} for each latent predictor.
#' @param latent_estimation_fn Optional function or list of functions used when
#'   \code{estimation_method = "custom"}.
#' @param ordinal Logical scalar or vector passed to \code{\link{lpmec_onerun}}.
#' @param mcmc_control List passed to \code{\link{lpmec_onerun}}.
#' @param min_split_correlation Minimum allowed split-half correlation. The
#'   correction is defined only for positive componentwise correlations.
#' @param conda_env Character string naming the conda environment for MCMC methods.
#' @param conda_env_required Logical indicating whether \code{conda_env} is required.
#'
#' @return A list of class \code{lpmec_multivariate_onerun} containing
#'   uncorrected OLS coefficients, uncorrected IV coefficients, corrected IV
#'   coefficients, split-half correlations, first-stage diagnostics, and latent
#'   score matrices.
#'
#' @export
lpmec_multivariate_onerun <- function(Y,
                                      observables,
                                      covariates = NULL,
                                      observables_groupings = NULL,
                                      make_observables_groupings = FALSE,
                                      estimation_method = "em",
                                      latent_estimation_fn = NULL,
                                      ordinal = FALSE,
                                      mcmc_control = list(
                                        backend = "pscl",
                                        n_samples_warmup = 500L,
                                        n_samples_mcmc = 1000L,
                                        batch_size = 512L,
                                        chain_method = "parallel",
                                        subsample_method = "full",
                                        anchor_parameter_id = NULL,
                                        n_thin_by = 1L,
                                        n_chains = 2L,
                                        outcome_prior = list(calibration = "data"),
                                        joint2_prior = list()),
                                      min_split_correlation = 0,
                                      conda_env = "lpmec",
                                      conda_env_required = FALSE) {
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
    stop("'Y' must contain only finite values for multivariate correction.")
  }
  if (!is.numeric(min_split_correlation) ||
      length(min_split_correlation) != 1L ||
      !is.finite(min_split_correlation) ||
      min_split_correlation < 0 ||
      min_split_correlation >= 1) {
    stop("'min_split_correlation' must be a single finite value in [0, 1).")
  }

  obs_list <- .lpmec_as_observables_list(observables)
  q <- length(obs_list)
  latent_names <- names(obs_list)
  if (is.null(latent_names) || any(!nzchar(latent_names))) {
    latent_names <- paste0("X", seq_len(q))
  }
  latent_names <- make.unique(latent_names)

  for (j in seq_len(q)) {
    if (!is.data.frame(obs_list[[j]]) && !is.matrix(obs_list[[j]])) {
      stop("'observables' element ", j, " must be a data.frame or matrix.")
    }
    if (nrow(obs_list[[j]]) != length(Y)) {
      stop("Number of rows in 'observables' element ", j,
           " must match length of 'Y'.")
    }
    if (ncol(obs_list[[j]]) < 4L) {
      stop("'observables' element ", j,
           " must have at least 4 columns to allow split-half estimation.")
    }
  }

  covariate_matrix <- .lpmec_prepare_covariates(covariates, length(Y))
  method_vec <- .lpmec_recycle_to_q(estimation_method, q, "estimation_method")
  make_group_vec <- .lpmec_recycle_to_q(
    make_observables_groupings,
    q,
    "make_observables_groupings"
  )
  ordinal_vec <- .lpmec_recycle_to_q(ordinal, q, "ordinal")
  grouping_list <- .lpmec_resolve_groupings(observables_groupings, obs_list)
  latent_fn_list <- .lpmec_resolve_latent_fns(latent_estimation_fn, q)

  x_est <- x_est1 <- x_est2 <- matrix(NA_real_, nrow = length(Y), ncol = q)
  colnames(x_est) <- colnames(x_est1) <- colnames(x_est2) <- latent_names
  scalar_runs <- vector("list", q)

  for (j in seq_len(q)) {
    scalar_runs[[j]] <- lpmec_onerun(
      Y = Y,
      observables = obs_list[[j]],
      observables_groupings = grouping_list[[j]],
      make_observables_groupings = make_group_vec[[j]],
      estimation_method = method_vec[[j]],
      latent_estimation_fn = latent_fn_list[[j]],
      ordinal = ordinal_vec[[j]],
      mcmc_control = mcmc_control,
      conda_env = conda_env,
      conda_env_required = conda_env_required
    )
    x_est[, j] <- as.numeric(scalar_runs[[j]]$x_est)
    x_est1[, j] <- as.numeric(scalar_runs[[j]]$x_est1)
    x_est2[, j] <- as.numeric(scalar_runs[[j]]$x_est2)
  }

  split_correlation <- vapply(seq_len(q), function(j) {
    stats::cor(x_est1[, j], x_est2[, j], use = "complete.obs")
  }, numeric(1L))
  names(split_correlation) <- latent_names
  if (any(!is.finite(split_correlation))) {
    stop("All split-half correlations must be finite.")
  }
  weak <- split_correlation <= min_split_correlation
  if (any(weak)) {
    stop(
      "Split-half correlations must be greater than 'min_split_correlation'. ",
      "Problematic latent predictor(s): ",
      paste(latent_names[weak], collapse = ", "),
      "."
    )
  }

  ols_fit <- .lpmec_fit_ols(Y, x_est, covariate_matrix)
  iv_fit_a <- .lpmec_fit_exact_iv(Y, x_est1, x_est2, covariate_matrix)
  iv_fit_b <- .lpmec_fit_exact_iv(Y, x_est2, x_est1, covariate_matrix)

  correction_factor <- sqrt(split_correlation)
  iv_coef_a <- iv_fit_a$latent_coef
  iv_coef_b <- iv_fit_b$latent_coef
  corrected_iv_coef_a <- correction_factor * iv_coef_a
  corrected_iv_coef_b <- correction_factor * iv_coef_b
  corrected_iv_coef <- (corrected_iv_coef_a + corrected_iv_coef_b) / 2
  iv_coef <- (iv_coef_a + iv_coef_b) / 2

  corrected_iv_coef_all <- c(
    corrected_iv_coef,
    (iv_fit_a$covariate_coef + iv_fit_b$covariate_coef) / 2
  )
  iv_coef_all <- c(iv_coef, (iv_fit_a$covariate_coef + iv_fit_b$covariate_coef) / 2)

  results <- list(
    n_latent = q,
    latent_names = latent_names,
    covariate_names = colnames(covariate_matrix),
    ols_coef = ols_fit$latent_coef,
    ols_se = ols_fit$latent_se,
    ols_coef_all = ols_fit$coef_all,
    ols_se_all = ols_fit$se_all,
    iv_coef_a = iv_coef_a,
    iv_coef_b = iv_coef_b,
    iv_coef = iv_coef,
    iv_coef_all = iv_coef_all,
    corrected_iv_coef_a = corrected_iv_coef_a,
    corrected_iv_coef_b = corrected_iv_coef_b,
    corrected_iv_coef = corrected_iv_coef,
    corrected_iv_coef_all = corrected_iv_coef_all,
    split_correlation = split_correlation,
    correction_factor = correction_factor,
    first_stage_fstat = .lpmec_first_stage_fstats(x_est1, x_est2, covariate_matrix),
    x_est = x_est,
    x_est1 = x_est1,
    x_est2 = x_est2,
    scalar_runs = scalar_runs
  )

  class(results) <- "lpmec_multivariate_onerun"
  results
}

#' Aggregated multivariate latent-predictor correction
#'
#' Runs \code{\link{lpmec_multivariate_onerun}} over repeated split-half
#' partitions and optional row bootstrap samples.
#'
#' @inheritParams lpmec_multivariate_onerun
#' @param n_boot Non-negative integer number of row-bootstrap iterations.
#' @param n_partition Positive integer number of split-half partitions per
#'   original or bootstrap sample.
#' @param partition_aggregation Aggregation strategy across partitions. See
#'   \code{\link{lpmec}}.
#' @param partition_aggregation_probs Quantile probabilities for winsorized or
#'   trimmed partition aggregation.
#' @param boot_basis Optional vector of indices or strata for row bootstrap.
#' @param return_intermediaries Logical. If \code{TRUE}, returns per-run
#'   coefficient matrices.
#'
#' @return A list of class \code{lpmec_multivariate} containing aggregated
#'   uncorrected OLS and corrected IV latent coefficients with bootstrap
#'   uncertainty summaries when \code{n_boot >= 1}.
#'
#' @export
lpmec_multivariate <- function(Y,
                               observables,
                               covariates = NULL,
                               observables_groupings = NULL,
                               make_observables_groupings = FALSE,
                               n_boot = 32L,
                               n_partition = 10L,
                               partition_aggregation = "median",
                               partition_aggregation_probs = c(0.01, 0.99),
                               boot_basis = seq_along(Y),
                               return_intermediaries = TRUE,
                               estimation_method = "em",
                               latent_estimation_fn = NULL,
                               ordinal = FALSE,
                               mcmc_control = list(
                                 backend = "pscl",
                                 n_samples_warmup = 500L,
                                 n_samples_mcmc = 1000L,
                                 batch_size = 512L,
                                 chain_method = "parallel",
                                 subsample_method = "full",
                                 anchor_parameter_id = NULL,
                                 n_thin_by = 1L,
                                 n_chains = 2L,
                                 outcome_prior = list(calibration = "data"),
                                 joint2_prior = list()),
                               min_split_correlation = 0,
                               conda_env = "lpmec",
                               conda_env_required = FALSE) {
  if (!is.numeric(n_boot) ||
      length(n_boot) != 1L ||
      !is.finite(n_boot) ||
      n_boot != floor(n_boot) ||
      n_boot < 0) {
    stop("'n_boot' must be a single non-negative integer.")
  }
  if (!is.numeric(n_partition) ||
      length(n_partition) != 1L ||
      !is.finite(n_partition) ||
      n_partition != floor(n_partition) ||
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

  obs_list <- .lpmec_as_observables_list(observables)
  the_sum_fxn <- .lpmec_resolve_partition_aggregation(
    partition_aggregation,
    partition_aggregation_probs
  )

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

      obs_boot <- lapply(obs_list, function(obs) obs[boot_indices, , drop = FALSE])
      cov_boot <- if (is.null(covariates)) {
        NULL
      } else {
        covariates[boot_indices, , drop = FALSE]
      }

      run <- lpmec_multivariate_onerun(
        Y = Y[boot_indices],
        observables = obs_boot,
        covariates = cov_boot,
        observables_groupings = observables_groupings,
        make_observables_groupings = make_observables_groupings,
        estimation_method = estimation_method,
        latent_estimation_fn = latent_estimation_fn,
        ordinal = ordinal,
        mcmc_control = mcmc_control,
        min_split_correlation = min_split_correlation,
        conda_env = conda_env,
        conda_env_required = conda_env_required
      )
      runs[[length(runs) + 1L]] <- run
      boot_ids <- c(boot_ids, boot_i)
      partition_ids <- c(partition_ids, partition_i)
    }
  }

  ols_runs <- .lpmec_runs_matrix(runs, "ols_coef")
  corrected_runs <- .lpmec_runs_matrix(runs, "corrected_iv_coef")
  iv_runs <- .lpmec_runs_matrix(runs, "iv_coef")
  ols_all_runs <- .lpmec_runs_matrix(runs, "ols_coef_all")
  corrected_all_runs <- .lpmec_runs_matrix(runs, "corrected_iv_coef_all")
  iv_all_runs <- .lpmec_runs_matrix(runs, "iv_coef_all")
  split_cor_runs <- .lpmec_runs_matrix(runs, "split_correlation")
  fstat_runs <- .lpmec_runs_matrix(runs, "first_stage_fstat")

  ols_by_boot <- .lpmec_aggregate_by_boot(ols_runs, boot_ids, the_sum_fxn)
  corrected_by_boot <- .lpmec_aggregate_by_boot(corrected_runs, boot_ids, the_sum_fxn)
  iv_by_boot <- .lpmec_aggregate_by_boot(iv_runs, boot_ids, the_sum_fxn)
  ols_all_by_boot <- .lpmec_aggregate_by_boot(ols_all_runs, boot_ids, the_sum_fxn)
  corrected_all_by_boot <- .lpmec_aggregate_by_boot(corrected_all_runs, boot_ids, the_sum_fxn)
  iv_all_by_boot <- .lpmec_aggregate_by_boot(iv_all_runs, boot_ids, the_sum_fxn)
  split_cor_by_boot <- .lpmec_aggregate_by_boot(split_cor_runs, boot_ids, the_sum_fxn)
  fstat_by_boot <- .lpmec_aggregate_by_boot(fstat_runs, boot_ids, the_sum_fxn)

  results <- list(
    n_latent = runs[[1L]]$n_latent,
    latent_names = runs[[1L]]$latent_names,
    covariate_names = runs[[1L]]$covariate_names,
    ols_coef = ols_by_boot[1L, ],
    ols_se = .lpmec_boot_sd(ols_by_boot, n_boot),
    ols_lower = .lpmec_boot_quantile(ols_by_boot, n_boot, 0.025),
    ols_upper = .lpmec_boot_quantile(ols_by_boot, n_boot, 0.975),
    iv_coef = iv_by_boot[1L, ],
    ols_coef_all = ols_all_by_boot[1L, ],
    ols_se_all = .lpmec_boot_sd(ols_all_by_boot, n_boot),
    iv_coef_all = iv_all_by_boot[1L, ],
    corrected_iv_coef = corrected_by_boot[1L, ],
    corrected_iv_se = .lpmec_boot_sd(corrected_by_boot, n_boot),
    corrected_iv_lower = .lpmec_boot_quantile(corrected_by_boot, n_boot, 0.025),
    corrected_iv_upper = .lpmec_boot_quantile(corrected_by_boot, n_boot, 0.975),
    corrected_iv_coef_all = corrected_all_by_boot[1L, ],
    corrected_iv_se_all = .lpmec_boot_sd(corrected_all_by_boot, n_boot),
    corrected_iv_lower_all = .lpmec_boot_quantile(corrected_all_by_boot, n_boot, 0.025),
    corrected_iv_upper_all = .lpmec_boot_quantile(corrected_all_by_boot, n_boot, 0.975),
    split_correlation = split_cor_by_boot[1L, ],
    first_stage_fstat = fstat_by_boot[1L, ],
    x_est = runs[[1L]]$x_est,
    x_est1 = runs[[1L]]$x_est1,
    x_est2 = runs[[1L]]$x_est2
  )

  if (return_intermediaries) {
    results <- c(results, list(
      Intermediary_BootIndex = boot_ids,
      Intermediary_PartitionIndex = partition_ids,
      Intermediary_ols_coef = ols_runs,
      Intermediary_iv_coef = iv_runs,
      Intermediary_corrected_iv_coef = corrected_runs,
      Intermediary_ols_coef_all = ols_all_runs,
      Intermediary_iv_coef_all = iv_all_runs,
      Intermediary_corrected_iv_coef_all = corrected_all_runs,
      Intermediary_split_correlation = split_cor_runs,
      Intermediary_first_stage_fstat = fstat_runs
    ))
  }

  class(results) <- "lpmec_multivariate"
  results
}

.lpmec_as_observables_list <- function(observables) {
  if (missing(observables) || is.null(observables)) {
    stop("'observables' is required and cannot be NULL.")
  }
  if (is.data.frame(observables) || is.matrix(observables)) {
    return(list(X1 = observables))
  }
  if (!is.list(observables) || length(observables) < 1L) {
    stop("'observables' must be a non-empty list, data.frame, or matrix.")
  }
  observables
}

.lpmec_recycle_to_q <- function(value, q, name) {
  if (length(value) == 1L) {
    return(rep(list(value), q))
  }
  if (length(value) != q) {
    stop("'", name, "' must have length 1 or length equal to the number of latent predictors.")
  }
  as.list(value)
}

.lpmec_resolve_groupings <- function(observables_groupings, obs_list) {
  q <- length(obs_list)
  if (is.null(observables_groupings)) {
    return(lapply(obs_list, colnames))
  }
  if (!is.list(observables_groupings) || is.data.frame(observables_groupings)) {
    if (q != 1L) {
      stop("'observables_groupings' must be a list when multiple latent predictors are supplied.")
    }
    return(list(observables_groupings))
  }
  if (length(observables_groupings) != q) {
    stop("'observables_groupings' must have one element per latent predictor.")
  }
  observables_groupings
}

.lpmec_resolve_latent_fns <- function(latent_estimation_fn, q) {
  if (is.null(latent_estimation_fn) || is.function(latent_estimation_fn)) {
    return(rep(list(latent_estimation_fn), q))
  }
  if (!is.list(latent_estimation_fn) || length(latent_estimation_fn) != q) {
    stop("'latent_estimation_fn' must be NULL, a function, or a list with one element per latent predictor.")
  }
  for (j in seq_len(q)) {
    if (!is.null(latent_estimation_fn[[j]]) && !is.function(latent_estimation_fn[[j]])) {
      stop("'latent_estimation_fn' element ", j, " must be NULL or a function.")
    }
  }
  latent_estimation_fn
}

.lpmec_prepare_covariates <- function(covariates, n) {
  if (is.null(covariates)) {
    out <- matrix(nrow = n, ncol = 0L)
    colnames(out) <- character(0L)
    return(out)
  }
  if (is.data.frame(covariates)) {
    if (nrow(covariates) != n) {
      stop("'covariates' must have the same number of rows as 'Y'.")
    }
    out <- stats::model.matrix(~ ., data = covariates)
    out <- out[, colnames(out) != "(Intercept)", drop = FALSE]
  } else if (is.matrix(covariates)) {
    if (nrow(covariates) != n) {
      stop("'covariates' must have the same number of rows as 'Y'.")
    }
    out <- covariates
  } else {
    stop("'covariates' must be NULL, a data.frame, or a matrix.")
  }
  out <- as.matrix(out)
  storage.mode(out) <- "double"
  if (any(!is.finite(out))) {
    stop("'covariates' must contain only finite numeric values after model-matrix expansion.")
  }
  if (is.null(colnames(out))) {
    colnames(out) <- paste0("C", seq_len(ncol(out)))
  }
  colnames(out) <- make.unique(colnames(out))
  out
}

.lpmec_design_matrix <- function(x, covariates) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  colnames(x) <- make.unique(colnames(x))
  out <- cbind("(Intercept)" = 1, x, covariates)
  if (any(!is.finite(out))) {
    stop("Latent score and covariate design matrices must contain only finite values.")
  }
  out
}

.lpmec_fit_ols <- function(Y, x, covariates) {
  design <- .lpmec_design_matrix(x, covariates)
  xtx <- crossprod(design)
  if (qr(xtx)$rank < ncol(xtx)) {
    stop("OLS design matrix is rank deficient.")
  }
  coef_all <- as.numeric(solve(xtx, crossprod(design, Y)))
  names(coef_all) <- colnames(design)
  resid <- as.numeric(Y - design %*% coef_all)
  df <- nrow(design) - ncol(design)
  if (df > 0L) {
    sigma2 <- sum(resid^2) / df
    se_all <- sqrt(diag(sigma2 * solve(xtx)))
  } else {
    se_all <- rep(NA_real_, length(coef_all))
  }
  names(se_all) <- names(coef_all)
  latent_idx <- seq.int(2L, ncol(x) + 1L)
  list(
    coef_all = coef_all[-1L],
    se_all = se_all[-1L],
    latent_coef = coef_all[latent_idx],
    latent_se = se_all[latent_idx]
  )
}

.lpmec_fit_exact_iv <- function(Y, x_regressor, x_instrument, covariates) {
  regressor <- .lpmec_design_matrix(x_regressor, covariates)
  instrument <- .lpmec_design_matrix(x_instrument, covariates)
  qzr <- crossprod(instrument, regressor) / nrow(regressor)
  qzy <- crossprod(instrument, Y) / nrow(regressor)
  if (qr(qzr)$rank < ncol(qzr)) {
    stop("IV moment matrix is rank deficient.")
  }
  theta <- as.numeric(solve(qzr, qzy))
  names(theta) <- colnames(regressor)
  latent_idx <- seq.int(2L, ncol(x_regressor) + 1L)
  covariate_idx <- if (ncol(covariates) > 0L) {
    seq.int(ncol(x_regressor) + 2L, length(theta))
  } else {
    integer(0L)
  }
  list(
    coef_all = theta[-1L],
    latent_coef = theta[latent_idx],
    covariate_coef = theta[covariate_idx]
  )
}

.lpmec_first_stage_fstats <- function(x_regressor, x_instrument, covariates) {
  out <- vapply(seq_len(ncol(x_regressor)), function(j) {
    design <- cbind("(Intercept)" = 1, instrument = x_instrument[, j], covariates)
    fit <- try(.lpmec_fit_ols(x_regressor[, j], design[, -1L, drop = FALSE], matrix(nrow = nrow(design), ncol = 0L)), silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(NA_real_)
    }
    se <- fit$latent_se[1L]
    coef <- fit$latent_coef[1L]
    if (!is.finite(se) || se <= 0) {
      return(NA_real_)
    }
    (coef / se)^2
  }, numeric(1L))
  names(out) <- colnames(x_regressor)
  out
}

.lpmec_runs_matrix <- function(runs, field) {
  out <- do.call(rbind, lapply(runs, function(run) run[[field]]))
  colnames(out) <- names(runs[[1L]][[field]])
  out
}

.lpmec_aggregate_by_boot <- function(values, boot_ids, aggregation_fn) {
  boots <- sort(unique(boot_ids))
  agg <- vapply(boots, function(boot) {
    apply(values[boot_ids == boot, , drop = FALSE], 2L, aggregation_fn)
  }, numeric(ncol(values)))
  # vapply() returns a plain vector when ncol(values) == 1L and a
  # ncol(values) x n_boots matrix otherwise; normalize to boots-by-columns.
  out <- if (is.matrix(agg)) t(agg) else matrix(agg, ncol = 1L)
  colnames(out) <- colnames(values)
  out
}

.lpmec_boot_sd <- function(values_by_boot, n_boot) {
  if (n_boot < 1L) {
    out <- rep(NA_real_, ncol(values_by_boot))
    names(out) <- colnames(values_by_boot)
    return(out)
  }
  out <- apply(values_by_boot[-1L, , drop = FALSE], 2L, stats::sd)
  names(out) <- colnames(values_by_boot)
  out
}

.lpmec_boot_quantile <- function(values_by_boot, n_boot, prob) {
  if (n_boot < 1L) {
    out <- rep(NA_real_, ncol(values_by_boot))
    names(out) <- colnames(values_by_boot)
    return(out)
  }
  out <- apply(values_by_boot[-1L, , drop = FALSE], 2L, function(x) {
    if (length(x) < 2L) {
      return(NA_real_)
    }
    stats::quantile(x, probs = prob, names = FALSE, na.rm = TRUE)
  })
  names(out) <- colnames(values_by_boot)
  out
}
