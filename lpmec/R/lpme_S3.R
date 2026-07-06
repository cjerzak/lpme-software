#' Summary method for lpmec_onerun objects
#'
#' Provides a summary of single-run LPMEC model results including OLS, IV,
#' and corrected coefficient estimates.
#'
#' @param object An object of class \code{lpmec_onerun} returned by \code{\link{lpmec_onerun}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame containing coefficient estimates and standard errors,
#'   returned invisibly. The data frame has rows for OLS, IV, Corrected IV,
#'   and Corrected OLS estimates.
#'
#' @seealso \code{\link{lpmec_onerun}}, \code{\link{print.lpmec_onerun}}, \code{\link{plot.lpmec_onerun}}
#'
#' @export
summary.lpmec_onerun <- function(object, ...) {
  coef_df <- data.frame(
    Estimate = c(object$ols_coef, object$iv_coef, object$corrected_iv_coef, object$corrected_ols_coef),
    SE = c(object$ols_se, object$iv_se, object$corrected_iv_se, object$corrected_ols_se),
    row.names = c("OLS", "IV", "Corrected IV", "Corrected OLS")
  )

  cat("Single-Run LPMEC Model Summary\n")
  cat("==============================\n")
  print(coef_df)
  invisible(coef_df)
}

#' Print method for lpmec_onerun objects
#'
#' Prints a concise summary of single-run LPMEC model results.
#'
#' @param x An object of class \code{lpmec_onerun} returned by \code{\link{lpmec_onerun}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object \code{x}, returned invisibly.
#'
#' @seealso \code{\link{lpmec_onerun}}, \code{\link{summary.lpmec_onerun}}, \code{\link{plot.lpmec_onerun}}
#'
#' @export
print.lpmec_onerun <- function(x, ...) {
  cat("Single-Run LPMEC Results\n")
  cat("------------------------\n")
  cat(sprintf("Uncorrected Coefficient (OLS): %.3f (SE: %.3f)\n", x$ols_coef, x$ols_se))
  cat(sprintf("Corrected Coefficient: %.3f (SE: %.3f)\n", x$corrected_ols_coef, x$corrected_ols_se))
  cat("Use summary() for detailed results.\n")
  invisible(x)
}

#' Plot method for lpmec_onerun objects
#'
#' Creates a scatter plot comparing the two split-half latent variable estimates.
#'
#' @param x An object of class \code{lpmec_onerun} returned by \code{\link{lpmec_onerun}}.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return No return value, called for side effects (creates a plot).
#'
#' @seealso \code{\link{lpmec_onerun}}, \code{\link{summary.lpmec_onerun}}, \code{\link{print.lpmec_onerun}}
#'
#' @export
plot.lpmec_onerun <- function(x, ...) {
  plot(x$x_est1, x$x_est2,
       xlab = "First Latent Estimate", ylab = "Second Latent Estimate",
       main = "Single-Run Latent Estimates", pch = 19, ...)
  abline(a = 0, b = 1, col = "blue", lty = 2)
}


#' Summary method for lpmec objects
#'
#' Provides a comprehensive summary of bootstrapped LPMEC model results including
#' OLS, IV, corrected, and Bayesian coefficient estimates with confidence intervals.
#'
#' @param object An object of class \code{lpmec} returned by \code{\link{lpmec}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame containing coefficient estimates, standard errors, and
#'   confidence intervals, returned invisibly. The data frame has rows for OLS,
#'   IV, Corrected IV, Corrected OLS, and Bayesian OLS estimates.
#'
#' @seealso \code{\link{lpmec}}, \code{\link{print.lpmec}}, \code{\link{plot.lpmec}}
#'
#' @export
summary.lpmec <- function(object, ...) {
  coef_df <- data.frame(
    Estimate = c(object$ols_coef, object$iv_coef, object$corrected_iv_coef,
                 object$corrected_ols_coef, object$bayesian_ols_coef_outer_normed,
                 object$bayesian_ols_coef_inner_normed),
    SE = c(object$ols_se, object$iv_se, object$corrected_iv_se,
           object$corrected_ols_se, object$bayesian_ols_se_outer_normed,
           object$bayesian_ols_se_inner_normed),
    CI_Lower = c(object$ols_lower, object$iv_lower, object$corrected_iv_lower,
                 object$corrected_ols_lower, object$bayesian_ols_lower_outer_normed,
                 object$bayesian_ols_lower_inner_normed),
    CI_Upper = c(object$ols_upper, object$iv_upper, object$corrected_iv_upper,
                 object$corrected_ols_upper, object$bayesian_ols_upper_outer_normed,
                 object$bayesian_ols_upper_inner_normed),
    row.names = c("OLS", "IV", "Corrected IV", "Corrected OLS",
                  "Bayesian OLS (Outer)", "Bayesian OLS (Inner)")
  )

  cat("Latent Predictor Measurement Error Correction (LPMEC) Model Summary\n")
  cat("====================================================================\n")
  if (!is.null(object$bootstrap_method)) {
    cat(sprintf(
      "Resampling: %s, m = %s (m/n = %.3f), CI = %s, replace = %s\n",
      object$bootstrap_method,
      object$boot_m,
      object$boot_m_ratio,
      object$boot_ci_type,
      object$boot_replace
    ))
    if (!is.null(object$bootstrap_success_rate) && is.finite(object$bootstrap_success_rate)) {
      cat(sprintf("Bootstrap success rate: %.3f\n", object$bootstrap_success_rate))
    }
  }
  print(coef_df)
  invisible(coef_df)
}

#' Print method for lpmec objects
#'
#' Prints a concise summary of bootstrapped LPMEC model results.
#'
#' @param x An object of class \code{lpmec} returned by \code{\link{lpmec}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object \code{x}, returned invisibly.
#'
#' @seealso \code{\link{lpmec}}, \code{\link{summary.lpmec}}, \code{\link{plot.lpmec}}
#'
#' @export
print.lpmec <- function(x, ...) {
  cat("Latent Predictor Measurement Error Correction (LPMEC) Model Results\n")
  cat("-------------------------------------------------------------------\n")
  if (!is.null(x$bootstrap_method)) {
    cat(sprintf(
      "Resampling: %s, m = %s, CI = %s\n",
      x$bootstrap_method,
      x$boot_m,
      x$boot_ci_type
    ))
  }
  cat(sprintf("Uncorrected Coefficient (OLS): %.3f (SE: %.3f)\n", x$ols_coef, x$ols_se))
  cat(sprintf("Corrected Coefficient: %.3f (SE: %.3f)\n", x$corrected_iv_coef, x$corrected_iv_se))
  cat(sprintf("Bayesian OLS (Outer): %.3f (SE: %.3f)\n", x$bayesian_ols_coef_outer_normed,
              x$bayesian_ols_se_outer_normed))
  cat("Use summary() for detailed results.\n")
  invisible(x)
}

#' Plot method for lpmec objects
#'
#' Creates visualizations of LPMEC model results. Can plot either the latent
#' variable estimates or the bootstrap distribution of coefficients.
#'
#' @param x An object of class \code{lpmec} returned by \code{\link{lpmec}}.
#' @param type Character string specifying the plot type. Either \code{"latent"}
#'   (default) for a scatter plot of split-half latent estimates, or
#'   \code{"coefficients"} for a density plot of bootstrap coefficient estimates.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}} or
#'   \code{\link[stats]{density}}.
#'
#' @return No return value, called for side effects (creates a plot).
#'
#' @seealso \code{\link{lpmec}}, \code{\link{summary.lpmec}}, \code{\link{print.lpmec}}
#'
#' @export
plot.lpmec <- function(x, type = "latent", ...) {
  if (type == "latent") {
    plot(x$x_est1, x$x_est2,
         xlab = "First Latent Estimate", ylab = "Second Latent Estimate",
         main = "Latent Variable Estimates", pch = 19, ...)
    abline(a = 0, b = 1, col = "red", lty = 2)
  } else if (type == "coefficients") {
    if (is.null(x$Intermediary_corrected_iv_coef) ||
        is.null(x$Intermediary_BootIndex)) {
      stop("Coefficient plots require lpmec(..., return_intermediaries = TRUE).")
    }
    boot_coefs <- as.numeric(x$Intermediary_corrected_iv_coef[
      x$Intermediary_BootIndex != 1
    ])
    boot_coefs <- boot_coefs[is.finite(boot_coefs)]
    if (length(boot_coefs) < 2L) {
      stop("Coefficient plots require bootstrap draws; run lpmec(..., n_boot >= 1, return_intermediaries = TRUE).")
    }
    plot(density(boot_coefs), main = "Bootstrap Distribution of Corrected IV Coefficient",
         xlab = "Coefficient Value", ...)
  } else {
    stop("Invalid plot type. Choose 'latent' or 'coefficients'.")
  }
}

#' Summary method for lpmec_multivariate_onerun objects
#'
#' @param object An object of class \code{lpmec_multivariate_onerun}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of latent-predictor coefficient estimates.
#'
#' @export
summary.lpmec_multivariate_onerun <- function(object, ...) {
  coef_df <- data.frame(
    OLS = object$ols_coef,
    IV = object$iv_coef,
    Corrected_IV = object$corrected_iv_coef,
    Split_Correlation = object$split_correlation,
    First_Stage_F = object$first_stage_fstat,
    row.names = object$latent_names
  )

  cat("Single-Run Multivariate LPMEC Summary\n")
  cat("=====================================\n")
  print(coef_df)
  invisible(coef_df)
}

#' Print method for lpmec_multivariate_onerun objects
#'
#' @param x An object of class \code{lpmec_multivariate_onerun}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object \code{x}, returned invisibly.
#'
#' @export
print.lpmec_multivariate_onerun <- function(x, ...) {
  cat("Single-Run Multivariate LPMEC Results\n")
  cat("-------------------------------------\n")
  cat(sprintf("Latent predictors: %s\n", paste(x$latent_names, collapse = ", ")))
  cat("Use summary() for coefficient details.\n")
  invisible(x)
}

#' Summary method for lpmec_multivariate objects
#'
#' @param object An object of class \code{lpmec_multivariate}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of aggregated latent-predictor coefficient estimates.
#'
#' @export
summary.lpmec_multivariate <- function(object, ...) {
  coef_df <- data.frame(
    OLS = object$ols_coef,
    OLS_SE = object$ols_se,
    Corrected_IV = object$corrected_iv_coef,
    Corrected_IV_SE = object$corrected_iv_se,
    Corrected_IV_Lower = object$corrected_iv_lower,
    Corrected_IV_Upper = object$corrected_iv_upper,
    Split_Correlation = object$split_correlation,
    First_Stage_F = object$first_stage_fstat,
    row.names = object$latent_names
  )

  cat("Multivariate LPMEC Summary\n")
  cat("==========================\n")
  print(coef_df)
  invisible(coef_df)
}

#' Print method for lpmec_multivariate objects
#'
#' @param x An object of class \code{lpmec_multivariate}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object \code{x}, returned invisibly.
#'
#' @export
print.lpmec_multivariate <- function(x, ...) {
  cat("Multivariate LPMEC Results\n")
  cat("--------------------------\n")
  cat(sprintf("Latent predictors: %s\n", paste(x$latent_names, collapse = ", ")))
  cat("Use summary() for coefficient details.\n")
  invisible(x)
}

#' Summary method for lpmec_panel_onerun objects
#'
#' @param object An object of class \code{lpmec_panel_onerun}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with one row per measure holding the naive OLS
#'   coefficient (with cluster-robust standard error), the headline corrected
#'   OLS coefficient with its [min, max] reliability-variant interval, the
#'   corrected IV coefficient, the pooled split and triad reliabilities, and
#'   the first-stage F statistic, returned invisibly.
#'
#' @export
summary.lpmec_panel_onerun <- function(object, ...) {
  coef_df <- data.frame(
    OLS = object$ols_coef,
    OLS_SE = object$ols_se,
    Corrected_OLS = object$corrected_ols_coef,
    Corrected_OLS_Lower = object$corrected_ols_lower,
    Corrected_OLS_Upper = object$corrected_ols_upper,
    Corrected_IV = object$corrected_iv_coef,
    Split_Correlation = object$split_correlation,
    Triad = object$triad,
    First_Stage_F = object$first_stage_fstat,
    row.names = object$measure_names
  )

  cat("Single-Run Panel LPMEC Summary\n")
  cat("==============================\n")
  cat(sprintf("Design: %s | Units: %s | Observations: %s\n",
              object$design, object$n_units, object$n_obs))
  print(coef_df)
  invisible(coef_df)
}

#' Print method for lpmec_panel_onerun objects
#'
#' @param x An object of class \code{lpmec_panel_onerun}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object \code{x}, returned invisibly.
#'
#' @export
print.lpmec_panel_onerun <- function(x, ...) {
  cat("Single-Run Panel LPMEC Results\n")
  cat("------------------------------\n")
  cat(sprintf("Design: %s | Units: %s | Observations: %s\n",
              x$design, x$n_units, x$n_obs))
  cat(sprintf("Measures: %s\n", paste(x$measure_names, collapse = ", ")))
  cat("Use summary() for coefficient details.\n")
  invisible(x)
}

#' Summary method for lpmec_panel objects
#'
#' @param object An object of class \code{lpmec_panel}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with one row per measure holding the naive OLS
#'   coefficient, the headline corrected OLS coefficient with its [min, max]
#'   reliability-variant interval, the corrected IV coefficient (with
#'   cluster-bootstrap standard errors), the pooled split and triad
#'   reliabilities, and the first-stage F statistic, returned invisibly.
#'
#' @export
summary.lpmec_panel <- function(object, ...) {
  coef_df <- data.frame(
    OLS = object$ols_coef,
    OLS_SE = object$ols_coef_se,
    Corrected_OLS = object$corrected_ols_coef,
    Corrected_OLS_SE = object$corrected_ols_coef_se,
    Corrected_OLS_Lower = object$corrected_ols_lower,
    Corrected_OLS_Upper = object$corrected_ols_upper,
    Corrected_IV = object$corrected_iv_coef,
    Corrected_IV_SE = object$corrected_iv_coef_se,
    Split_Correlation = object$split_correlation,
    Triad = object$triad,
    First_Stage_F = object$first_stage_fstat,
    row.names = object$measure_names
  )

  cat("Panel LPMEC Summary\n")
  cat("===================\n")
  cat(sprintf(
    "Design: %s | Units: %s | Observations: %s | Bootstrap replications: %s\n",
    object$design, object$n_units, object$n_obs, object$n_boot
  ))
  print(coef_df)
  invisible(coef_df)
}

#' Print method for lpmec_panel objects
#'
#' @param x An object of class \code{lpmec_panel}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object \code{x}, returned invisibly.
#'
#' @export
print.lpmec_panel <- function(x, ...) {
  cat("Panel LPMEC Results\n")
  cat("-------------------\n")
  cat(sprintf("Design: %s | Units: %s | Observations: %s\n",
              x$design, x$n_units, x$n_obs))
  cat(sprintf("Measures: %s\n", paste(x$measure_names, collapse = ", ")))
  cat(sprintf("Bootstrap replications: %s | Partitions: %s\n",
              x$n_boot, x$n_partition))
  cat("Use summary() for coefficient details.\n")
  invisible(x)
}

#' Plot method for lpmec_panel objects
#'
#' Creates a scatter plot of the two pooled half scores of the first measure
#' with split halves; when no measure has half scores, falls back to a
#' scatter plot of the first two measure scores.
#'
#' @param x An object of class \code{lpmec_panel}.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return No return value, called for side effects (creates a plot).
#'
#' @export
plot.lpmec_panel <- function(x, ...) {
  split_available <- vapply(x$measure_names, function(m) {
    sum(is.finite(x$x_est1[, m]) & is.finite(x$x_est2[, m])) >= 2L
  }, logical(1L))
  if (any(split_available)) {
    measure <- x$measure_names[split_available][1L]
    plot(x$x_est1[, measure], x$x_est2[, measure],
         xlab = "First Half Score", ylab = "Second Half Score",
         main = sprintf("Panel Measure Half Scores (%s)", measure),
         pch = 19, ...)
  } else if (x$n_measures >= 2L) {
    plot(x$x_est[, 1L], x$x_est[, 2L],
         xlab = sprintf("Measure Score (%s)", x$measure_names[1L]),
         ylab = sprintf("Measure Score (%s)", x$measure_names[2L]),
         main = "Panel Measure Scores", pch = 19, ...)
  } else {
    stop("Plotting requires a measure with half scores or at least two measures.")
  }
  abline(a = 0, b = 1, col = "blue", lty = 2)
}

#' Summary method for lpmec_moderator_onerun objects
#'
#' @param object An object of class \code{lpmec_moderator_onerun}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of naive and corrected interaction-regression
#'   coefficient estimates, returned invisibly. When three or more measures
#'   are supplied, the per-measure triad corrections are also printed.
#'
#' @export
summary.lpmec_moderator_onerun <- function(object, ...) {
  coef_df <- data.frame(
    Estimate = c(object$treatment_coef, object$main_coef,
                 object$interaction_coef, object$corrected_main_coef,
                 object$corrected_interaction_coef),
    SE = c(object$treatment_se, object$main_se, object$interaction_se,
           NA_real_, NA_real_),
    row.names = c("Treatment", "Moderator", "Interaction",
                  "Corrected Moderator", "Corrected Interaction")
  )

  cat("Single-Run Latent-Moderator LPMEC Summary\n")
  cat("=========================================\n")
  cat(sprintf(
    "Reliability: rho_half = %.3f, SB factor = %s, rho_score = %.3f\n",
    object$rho_half, object$sb_factor, object$rho_score
  ))
  print(coef_df)
  if (!is.null(object$per_measure)) {
    cat("\nPer-measure triad corrections:\n")
    print(object$per_measure)
  }
  invisible(coef_df)
}

#' Print method for lpmec_moderator_onerun objects
#'
#' @param x An object of class \code{lpmec_moderator_onerun}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object \code{x}, returned invisibly.
#'
#' @export
print.lpmec_moderator_onerun <- function(x, ...) {
  cat("Single-Run Latent-Moderator LPMEC Results\n")
  cat("-----------------------------------------\n")
  cat(sprintf("Naive Interaction: %.3f (SE: %.3f)\n",
              x$interaction_coef, x$interaction_se))
  cat(sprintf("Corrected Interaction: %.3f (rho_score: %.3f)\n",
              x$corrected_interaction_coef, x$rho_score))
  cat("Use summary() for detailed results.\n")
  invisible(x)
}

#' Summary method for lpmec_moderator objects
#'
#' @param object An object of class \code{lpmec_moderator}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of naive and corrected interaction-regression
#'   coefficient estimates and reliability quantities with bootstrap
#'   standard errors and percentile confidence intervals, returned
#'   invisibly.
#'
#' @export
summary.lpmec_moderator <- function(object, ...) {
  coef_df <- data.frame(
    Estimate = c(object$treatment_coef, object$main_coef,
                 object$interaction_coef, object$corrected_main_coef,
                 object$corrected_interaction_coef, object$rho_half,
                 object$rho_score),
    SE = c(object$treatment_se, object$main_se, object$interaction_se,
           object$corrected_main_se, object$corrected_interaction_se,
           object$rho_half_se, object$rho_score_se),
    CI_Lower = c(object$treatment_lower, object$main_lower,
                 object$interaction_lower, object$corrected_main_lower,
                 object$corrected_interaction_lower, object$rho_half_lower,
                 object$rho_score_lower),
    CI_Upper = c(object$treatment_upper, object$main_upper,
                 object$interaction_upper, object$corrected_main_upper,
                 object$corrected_interaction_upper, object$rho_half_upper,
                 object$rho_score_upper),
    row.names = c("Treatment", "Moderator", "Interaction",
                  "Corrected Moderator", "Corrected Interaction",
                  "rho_half", "rho_score")
  )

  cat("Latent-Moderator LPMEC Summary\n")
  cat("==============================\n")
  cat(sprintf(
    "Resampling: %s bootstrap replication(s), %s partition(s), SB factor = %s\n",
    object$n_boot, object$n_partition, object$sb_factor
  ))
  print(coef_df)
  invisible(coef_df)
}

#' Print method for lpmec_moderator objects
#'
#' @param x An object of class \code{lpmec_moderator}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object \code{x}, returned invisibly.
#'
#' @export
print.lpmec_moderator <- function(x, ...) {
  cat("Latent-Moderator LPMEC Results\n")
  cat("------------------------------\n")
  cat(sprintf("Naive Interaction: %.3f (SE: %.3f)\n",
              x$interaction_coef, x$interaction_se))
  cat(sprintf("Corrected Interaction: %.3f (SE: %.3f)\n",
              x$corrected_interaction_coef, x$corrected_interaction_se))
  cat(sprintf("Score Reliability (rho_score): %.3f (SE: %.3f)\n",
              x$rho_score, x$rho_score_se))
  cat("Use summary() for detailed results.\n")
  invisible(x)
}

#' Plot method for lpmec_moderator objects
#'
#' Creates a scatter plot of the two pooled half scores of the first measure
#' with split halves; when no measure has half scores, falls back to a
#' scatter plot of the first two measure scores.
#'
#' @param x An object of class \code{lpmec_moderator}.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return No return value, called for side effects (creates a plot).
#'
#' @export
plot.lpmec_moderator <- function(x, ...) {
  split_available <- vapply(x$measure_names, function(m) {
    sum(is.finite(x$x_est1[, m]) & is.finite(x$x_est2[, m])) >= 2L
  }, logical(1L))
  if (any(split_available)) {
    measure <- x$measure_names[split_available][1L]
    plot(x$x_est1[, measure], x$x_est2[, measure],
         xlab = "First Half Score", ylab = "Second Half Score",
         main = sprintf("Moderator Half Scores (%s)", measure),
         pch = 19, ...)
  } else if (x$n_measures >= 2L) {
    plot(x$x_est[, 1L], x$x_est[, 2L],
         xlab = sprintf("Measure Score (%s)", x$measure_names[1L]),
         ylab = sprintf("Measure Score (%s)", x$measure_names[2L]),
         main = "Moderator Measure Scores", pch = 19, ...)
  } else {
    stop("Plotting requires a measure with half scores or at least two measures.")
  }
  abline(a = 0, b = 1, col = "blue", lty = 2)
}

#' Print method for lpmec_reliability_bounds objects
#'
#' @param x An object of class \code{lpmec_reliability_bounds}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object \code{x}, returned invisibly.
#'
#' @export
print.lpmec_reliability_bounds <- function(x, ...) {
  cat("LPMEC Reliability Bounds\n")
  cat("------------------------\n")
  cat(sprintf("Measures: %s\n", paste(x$measure_names, collapse = ", ")))
  cat(sprintf("Designs: %s\n", paste(x$designs, collapse = ", ")))
  if (!is.null(x$n_boot) && x$n_boot > 0L) {
    cat(sprintf("Bootstrap replications: %s (failed: %s)\n",
                x$n_boot, x$n_boot_failed))
  }
  reliability <- x$reliability
  numeric_columns <- vapply(reliability, is.numeric, logical(1L))
  reliability[numeric_columns] <- lapply(
    reliability[numeric_columns], round, digits = 3
  )
  print(reliability)
  invisible(x)
}
