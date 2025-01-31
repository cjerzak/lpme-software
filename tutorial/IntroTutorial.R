{
  # install.packages("~/Documents/lpme-software/lpme", repos = NULL, type = "source", force = FALSE)
  # Install via GitHub
  # devtools::install_github("cjerzak/lpme-software/lpme")
  
  options(error = NULL)
  n <- 5000  # Number of observations
  d <- 24    # Number of observable indicators
  
  # 1) Generate latent ability
  x_true <- rnorm(n)  # 
  
  # 2) Generate item parameters:
  #    - All discriminations are positive to ensure positive correlations
  #    - Difficulties are varied
  discrimination <- runif(d, 0.5, 2)  # a_j
  difficulty     <- runif(d, -1, 1)   # b_j
  
  # 3) For each item j, compute the probability of responding "1" via the 2PL model
  #      P(X_{ij} = 1) = 1 / [1 + exp(-a_j (theta_i - b_j))]
  #    We use plogis() which is 1 / (1 + exp(-z)).
  P <- sapply(1:d, function(j) {
    plogis(discrimination[j] * (x_true - difficulty[j]))
  })
  
  # 4) Generate binary responses for each item
  ObservablesMat <- sapply(1:d, function(j) {
    rbinom(n, 1, P[, j])
  })
  
  # 5) Generate a continuous outcome correlated with the latent trait
  Yobs <- 0.4 * x_true + rnorm(n, sd = 0.35)
  
  library( lpme )
  
  # Run bootstrapped analysis
  results <- lpme(
    Y = Yobs,
    observables = as.data.frame(ObservablesMat),
    n_boot = 6L,      # Reduced for demonstration
    n_partition = 5L,  # Reduced for demonstration
    estimation_method = "emIRT"
  )
  
  # Compare estimates
  print(results)
  summary(results)
  
  # Visualization
  plot(results)
  
  # Bayesian MCMC example (commented out)
  # lpme::build_backend() # build backend if needed
  t0 <- Sys.time()
  mcmc_results <- lpme(
      Y = Yobs,
      observables = as.data.frame(ObservablesMat),
      n_boot = 0L,      # Reduced for demonstration
      n_partition = 2L, # Reduced for demonstration
      estimation_method = "MCMC",
      mcmc_control = list(
                #backend = "numpyro",  
                backend = "pscl",  
                n_samples_warmup = 500L, n_samples_mcmc = 1000L, subsample_method = "full", 
                #n_samples_warmup = 500L, n_samples_mcmc = 1000L, subsample_method = "batch", batch_size = 128L, 
                chain_method = "sequential", 
                n_thin_by = 1L, 
                n_chains = 2L), 
      conda_env = "lpme"  # Specify your conda environment
      #conda_env = "jax_gpu"  # Specify your conda environment
    )
  t1 <- Sys.time()
  print(sprintf("MCMC runtime of %.4f mins", difftime(t1,t0, units = "mins")))
  print(mcmc_results)
  summary(mcmc_results)
  
  # Visualization
  plot(mcmc_results)
  
  # compare 
  summary(results)
  summary(mcmc_results)
  
  
}
