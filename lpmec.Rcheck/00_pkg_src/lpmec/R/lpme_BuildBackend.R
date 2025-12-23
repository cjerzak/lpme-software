#' A function to build the environment for lpmec. Builds a conda environment in which
#' 'JAX', 'numpyro', and 'numpy' are installed. Users can also create a conda
#' environment where 'JAX' and 'numpy' are installed themselves.
#'
#' @param conda_env (default = `"lpmec"`) Name of the conda environment in which to place the backends.
#' @param conda (default = `auto`) The path to a conda executable. Using `"auto"` allows reticulate to attempt to automatically find an appropriate conda binary.

#' @return Invisibly returns NULL; this function is used for its side effects
#' of creating and configuring a conda environment for `lpmec`.
#' This function requires an Internet connection.
#' You can find out a list of conda Python paths via: `Sys.which("python")`
#'
#' @examples
#' \dontrun{
#' # Create a conda environment named "lpmec"
#' # and install the required Python packages (jax, numpy, etc.)
#' build_backend(conda_env = "lpmec", conda = "auto")
#'
#' # If you want to specify a particular conda path:
#' # build_backend(conda_env = "lpmec", conda = "/usr/local/bin/conda")
#' }
#'
#' @export
#' @md

build_backend <- function(conda_env = "lpmec", conda = "auto"){
  # Create a new conda environment
  reticulate::conda_create(envname = conda_env,
                           conda = conda,
                           python_version = "3.13")
  
  # Define packages to install 
  Packages2Install <- c("jax", 
                        "numpy",
                        "numpyro")
  
  # Install packages 
  reticulate::py_install(Packages2Install, conda = conda, pip = TRUE, envname = conda_env)
}
