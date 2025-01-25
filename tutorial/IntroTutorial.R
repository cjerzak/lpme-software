{
  # install.packages("~/Documents/lpme-software/lpme", repos = NULL, type = "source", force = FALSE)
  # Install from GitHub
  # devtools::install_github("cjerzak/lpme-software/lpme")
  
  set.seed(123); options( error = NULL )
  n <- 500  # Number of observations
  d <- 10    # Number of observable indicators
  
  # Generate latent variable and observed outcomes
  x_true <- rnorm(n)
  Yobs <- 0.4 * x_true + rnorm(n, sd = 0.35)
  
  # Generate binary indicators of latent variable
  ObservablesMat <- sapply(1:d, function(j) {
    p <- pnorm(0.5 * x_true + rnorm(n, sd = 0.5))
    rbinom(n, 1, p)
  })
  
  library(lpme)
  
  # Run bootstrapped analysis
  results <- lpme(
    Y = Yobs,
    observables = as.data.frame(ObservablesMat),
    n_boot = 50,      # Reduced for demonstration
    n_partition = 5,  # Reduced for demonstration
    estimation_method = "emIRT"
  )
  
  # Compare estimates
  print(results)
  summary(results)
  
  # Visualization
  plot(results)
  
  # Bayesian MCMC example (commented out)
  # lpme::build_backend() # build backend if needed
  mcmc_results <- lpme(
      Y = Yobs,
      observables = as.data.frame(ObservablesMat),
      n_boot = 0L, 
      n_partition = 3L, 
      estimation_method = "MCMC",
      conda_env = "lpme"  # Specify your conda environment
    )
}
