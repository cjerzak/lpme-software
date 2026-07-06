make_lpmec_test_data <- function(n = 80L, p = 6L, seed = 123L) {
  set.seed(seed)
  list(
    Y = rnorm(n),
    obs = as.data.frame(matrix(sample(c(0, 1), n * p, replace = TRUE), ncol = p))
  )
}

# Balanced panel DGP (port of V2/Code/sim_paper/00_common.R make_panel):
# latent X_it = A_i + W_it with Var(A) = 1 - s and W_it stationary AR(1)(phi)
# with marginal variance s, so pooled Var(X) = 1; split scores t1/t2 = X + U
# with independent errors of variance sig2U each; Y = beta * X + unit FE +
# time FE + noise calibrated to hit the target within R^2.
make_panel_test_data <- function(n_units = 30L,
                                 n_periods = 12L,
                                 s = 0.5,
                                 sig2U = 0.5,
                                 beta = 0.4,
                                 phi = 0.95,
                                 r2_within = 0.15,
                                 sig2_alpha = 1,
                                 sig2_delta = 0.25,
                                 seed = 42L) {
  set.seed(seed)
  ar1_matrix <- function(N, T_periods, phi, v) {
    M <- matrix(0, N, T_periods)
    M[, 1] <- rnorm(N, 0, sqrt(v))
    if (T_periods > 1) {
      innovation_sd <- sqrt(v * (1 - phi^2))
      for (t in 2:T_periods) {
        M[, t] <- phi * M[, t - 1] + rnorm(N, 0, innovation_sd)
      }
    }
    M
  }
  ar1_within_var <- function(phi, v, T_periods) {
    g <- v * phi^abs(outer(seq_len(T_periods), seq_len(T_periods), "-"))
    v - mean(g)
  }
  sig2_eps <- beta^2 * ar1_within_var(phi, s, n_periods) *
    (1 - r2_within) / r2_within

  A <- rnorm(n_units, 0, sqrt(1 - s))
  X <- A + ar1_matrix(n_units, n_periods, phi, s)  # A recycles over columns
  alpha_unit <- rnorm(n_units, 0, sqrt(sig2_alpha))
  delta_time <- matrix(rnorm(n_periods, 0, sqrt(sig2_delta)),
                       n_units, n_periods, byrow = TRUE)
  eps <- matrix(rnorm(n_units * n_periods, 0, sqrt(sig2_eps)),
                n_units, n_periods)
  Y <- beta * X + alpha_unit + delta_time + eps
  t1 <- X + matrix(rnorm(n_units * n_periods, 0, sqrt(sig2U)),
                   n_units, n_periods)
  t2 <- X + matrix(rnorm(n_units * n_periods, 0, sqrt(sig2U)),
                   n_units, n_periods)

  list(
    Y = as.vector(Y),
    unit = rep(seq_len(n_units), times = n_periods),
    time = rep(seq_len(n_periods), each = n_units),
    X = as.vector(X),
    t1 = as.vector(t1),
    t2 = as.vector(t2),
    Y_mat = Y,
    X_mat = X,
    t1_mat = t1,
    t2_mat = t2,
    n_units = n_units,
    n_periods = n_periods,
    beta = beta,
    phi = phi,
    s = s,
    sig2U = sig2U,
    r2_within = r2_within
  )
}

# Latent-moderator experiment DGP with a 2PL-lite probit item battery (port of
# the calibrated simulator in V2/Code/moderator/00_helpers.R):
# X ~ N(0, 1); treatment ~ Bernoulli(0.5);
# Y = b_treatment * T + b_moderator * X + b_interaction * T * X + N(0, 1);
# P(item_j = 1 | X) = pnorm(difficulty_j + discrimination_j * X).
make_moderator_test_data <- function(n = 400L,
                                     n_items = 4L,
                                     b_treatment = 0.2,
                                     b_moderator = 0.2,
                                     b_interaction = 0.3,
                                     discrimination = NULL,
                                     difficulty = NULL,
                                     seed = 202L) {
  set.seed(seed)
  if (is.null(discrimination)) {
    discrimination <- rep(c(1.2, 1.0, 0.8, 1.1), length.out = n_items)
  }
  if (is.null(difficulty)) {
    difficulty <- rep(c(-0.4, 0.0, 0.3, -0.1), length.out = n_items)
  }
  X <- rnorm(n)
  treatment <- rbinom(n, 1, 0.5)
  Y <- b_treatment * treatment + b_moderator * X +
    b_interaction * treatment * X + rnorm(n)
  item_probs <- stats::pnorm(
    sweep(outer(X, discrimination), 2, difficulty, "+")
  )
  items <- matrix(rbinom(n * n_items, 1, as.vector(item_probs)), n, n_items)
  colnames(items) <- paste0("item", seq_len(n_items))

  list(
    Y = Y,
    treatment = treatment,
    X = X,
    items = items,
    b_treatment = b_treatment,
    b_moderator = b_moderator,
    b_interaction = b_interaction,
    discrimination = discrimination,
    difficulty = difficulty
  )
}
