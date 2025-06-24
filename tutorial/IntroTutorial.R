{
  # Local install for development team 
  # install.packages("~/Documents/lpme-software/lpme", repos = NULL, type = "source", force = FALSE)
  
  # Install via GitHub
  # devtools::install_github("cjerzak/lpme-software/lpme")
  
  options(error = NULL)
  n <- 600  # Number of observations
  d <- 6    # Number of observable indicators
  
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
    estimation_method = "em"
  )
  
  # Compare estimates
  print(results)
  summary(results)
  
  # Visualization
  plot(results)
  
  # Advanced use
  if(TRUE == FALSE){ 
  # Try out the Bayesian methods 
  # lpme::build_backend() # build backend if needed
  BayesBackend <- "pscl"
  # BayesBackend <- "numpyro"
  mcmc_results <- lpme(
      Y = Yobs,
      observables = as.data.frame(ObservablesMat),
      n_boot = 0L,      # Reduced for demonstration
      n_partition = 2L, # Reduced for demonstration
      estimation_method = "mcmc",
      mcmc_control = list(
                backend = BayesBackend,  
                n_samples_warmup = 1000L, n_samples_mcmc = 1000L, subsample_method = "full", 
                chain_method = "sequential", 
                n_thin_by = 1L, 
                n_chains = 2L), 
      conda_env = "lpme"  # Specify your conda environment, used in this condition, backend="numpyro"
  )
  mcmc_overimputation_results <- lpme(
    Y = Yobs,
    observables = as.data.frame(ObservablesMat),
    n_boot = 0L,      # Reduced for demonstration
    n_partition = 2L, # Reduced for demonstration
    estimation_method = "mcmc_overimputation",
    mcmc_control = list(
      backend = BayesBackend,   
      n_samples_warmup = 1000L, n_samples_mcmc = 1000L, subsample_method = "full", 
      chain_method = "sequential", 
      n_thin_by = 1L, 
      n_chains = 2L), 
  conda_env = "lpme"  # Specify your conda environment, used in this condition, backend="numpyro"
  )
  mcmc_joint_results <- lpme(
    Y = Yobs,
    observables = as.data.frame(ObservablesMat),
    n_boot = 0L,      # Reduced for demonstration
    n_partition = 1L, # Reduced for demonstration
    estimation_method = "mcmc_joint",
    mcmc_control = list(
      backend = BayesBackend,  
      n_samples_warmup = 1000L, n_samples_mcmc = 1000L, subsample_method = "full", 
      chain_method = "sequential", 
      n_thin_by = 1L, 
      n_chains = 2L), 
    conda_env = "lpme"  # Specify your conda environment, used in this condition, backend="numpyro"
  )

  # compare 
  summary(results)
  summary(mcmc_results)
  summary(mcmc_joint_results)
  summary(mcmc_overimputation_results)
  }

}
