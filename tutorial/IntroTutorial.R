{
  # install.packages("~/Documents/lpme-software/lpme", repos = NULL, type = "source", force = FALSE)
  # Install from GitHub
  # devtools::install_github("cjerzak/lpme-software/lpme")
  
  set.seed(123); options( error = NULL )
  n <- 500  # Number of observations
  d <- 8    # Number of observable indicators
  
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
    n_boot = 8L,      # Reduced for demonstration
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
  mcmc_results <- lpme(
      Y = Yobs,
      observables = as.data.frame(ObservablesMat),
      n_boot = 0L,      # Reduced for demonstration
      n_partition = 3L, # Reduced for demonstration
      estimation_method = "MCMC",
      conda_env = "lpme"  # Specify your conda environment
    )
}

