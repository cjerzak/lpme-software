#' Infer orientation signs for each observable indicator
#'
#' This helper analyzes observable indicators and returns a numeric vector
#' of \code{1} or \code{-1} for use with the \code{orientation_signs}
#' argument in \code{lpme}. Each sign is chosen so that the correlation
#' between the oriented indicator and either the outcome \code{Y} or the
#' first principal component of the indicators is positive.
#' @importFrom stats cor prcomp
#'
#' @param Y Numeric outcome vector. Only used when \code{method = "Y"}.
#' @param observables A matrix or data frame of binary observable indicators.
#' @param method Character string specifying how to orient the indicators.
#'   \describe{
#'     \item{\code{"Y"}}{orient each indicator so that its correlation with
#'       \code{Y} is positive.}
#'     \item{\code{"PC1"}}{orient each indicator so that its correlation with
#'       the first principal component of \code{observables} is positive.}
#'   }
#'   Default is \code{"Y"}.
#'
#' @return A numeric vector of length \code{ncol(observables)} containing
#'   \code{1} or \code{-1}.
#' @examples
#' set.seed(1)
#' Y <- rnorm(10)
#' obs <- data.frame(matrix(sample(c(0,1), 20, replace = TRUE), ncol = 2))
#' infer_orientation_signs(Y, obs)
#' @export
infer_orientation_signs <- function(Y, observables, method = c("Y", "PC1")) {
  method <- match.arg(method)
  # allow for missing values when validating the observables
  obs_vals <- na.omit(unlist(observables))
  if (!all(obs_vals %in% c(0, 1))) {
    stop("infer_orientation_signs currently supports only binary observables")
  }
  if (method == "PC1") {
    pc1 <- prcomp(observables, scale. = TRUE)$x[, 1]
    cor_target <- pc1
  } else {
    if (missing(Y)) {
      stop("Y must be provided when method = 'Y'")
    }
    cor_target <- Y
  }
  signs <- apply(observables, 2, function(x) {
    cval <- suppressWarnings(cor(x, cor_target, use = "pairwise.complete.obs"))
    ifelse(is.na(cval) || cval >= 0, 1, -1)
  })
  names(signs) <- colnames(observables)
  as.numeric(signs)
}
