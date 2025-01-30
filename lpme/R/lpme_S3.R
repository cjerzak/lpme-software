#' @export
summary.lpme_onerun <- function(object, ...) {
  coef_df <- data.frame(
    Estimate = c(object$ols_coef, object$iv_coef, object$corrected_iv_coef, object$corrected_ols_coef),
    SE = c(object$ols_se, object$iv_se, object$corrected_iv_se, object$corrected_ols_se),
    row.names = c("OLS", "IV", "Corrected IV", "Corrected OLS")
  )
  
  cat("Single-Run LPME Model Summary\n")
  cat("============================\n")
  print(coef_df)
  invisible(coef_df)
}

#' @export
print.lpme_onerun <- function(x, ...) {
  cat("Single-Run LPME Results\n")
  cat("-----------------------\n")
  cat(sprintf("Uncorrected Coefficient (OLS): %.3f (SE: %.3f)\n", x$ols_coef, x$ols_se))
  cat(sprintf("Corrected Coefficient: %.3f (SE: %.3f)\n", x$ols_coef, x$ols_se))
  cat("Use summary() for detailed results.\n")
}

#' @export
plot.lpme_onerun <- function(x, ...) {
  plot(x$x_est1, x$x_est2,
       xlab = "First Latent Estimate", ylab = "Second Latent Estimate",
       main = "Single-Run Latent Estimates", pch = 19, ...)
  abline(a = 0, b = 1, col = "blue", lty = 2)
}


#' @export
summary.lpme <- function(object, ...) {
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
  
  cat("Latent Profile Measurement Error (LPME) Model Summary\n")
  cat("====================================================\n")
  print(coef_df)
  invisible(coef_df)
}

#' @export
print.lpme <- function(x, ...) {
  cat("Latent Profile Measurement Error (LPME) Model Results\n")
  cat("----------------------------------------------------\n")
  cat(sprintf("Uncorrected Coefficient (OLS): %.3f (SE: %.3f)\n", x$ols_coef, x$ols_se))
  cat(sprintf("Corrected Coefficient: %.3f (SE: %.3f)\n", x$corrected_iv_coef, x$corrected_iv_se))
  cat(sprintf("Bayesian OLS (Outer): %.3f (SE: %.3f)\n", x$bayesian_ols_coef_outer_normed, 
              x$bayesian_ols_se_outer_normed))
  cat("Use summary() for detailed results.\n")
}

#' @export
plot.lpme <- function(x, type = "latent", ...) {
  if (type == "latent") {
    plot(x$x_est1, x$x_est2, 
         xlab = "First Latent Estimate", ylab = "Second Latent Estimate",
         main = "Latent Variable Estimates", pch = 19, ...)
    abline(a = 0, b = 1, col = "red", lty = 2)
  } else if (type == "coefficients") {
    boot_coefs <- x$Intermediary_corrected_iv_coef[, -1]
    plot(density(boot_coefs), main = "Bootstrap Distribution of Corrected IV Coefficient",
         xlab = "Coefficient Value", ...)
  } else {
    stop("Invalid plot type. Choose 'latent' or 'coefficients'.")
  }
}