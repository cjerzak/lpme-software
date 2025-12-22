#' Initialize JAX/numpyro Python Environment
#'
#' Internal function to initialize the JAX and numpyro Python packages
#' for MCMC estimation via reticulate.
#'
#' @param conda_env Name of conda environment to use
#' @param conda_env_required Whether the conda environment is required
#'
#' @noRd
#' @keywords internal
initialize_jax <- function(conda_env = "lpme",
                           conda_env_required = TRUE) {

  # Load reticulate (Declared in Imports: in DESCRIPTION)
  reticulate::use_condaenv(condaenv = conda_env, required = conda_env_required)
  
  # Import Python packages once, storing them in lpme_env
  if (!exists("jax", envir = lpme_env, inherits = FALSE)) {
    lpme_env$jax <- reticulate::import("jax")
    lpme_env$jnp <- reticulate::import("jax.numpy")
    lpme_env$np  <- reticulate::import("numpy")
    lpme_env$random  <- reticulate::import("jax.random")
    lpme_env$numpyro  <- reticulate::import("numpyro")
    lpme_env$dist  <- reticulate::import("numpyro.distributions")
  }
} 
