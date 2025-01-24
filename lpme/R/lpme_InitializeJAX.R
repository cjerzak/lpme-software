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
