## ----setup, eval=FALSE--------------------------------------------------------
# # Install lpmec from source (replace with appropriate installation method)
# # devtools::install_github("cjerzak/lpmec-software/lpmec")

## ----eval=TRUE----------------------------------------------------------------
set.seed(123)
n <- 1000  # Number of observations
d <- 10    # Number of observable indicators

# Generate latent variable and observed outcomes
x_true <- rnorm(n)
Yobs <- 0.4 * x_true + rnorm(n, sd = 0.35)

# Generate binary indicators of latent variable
ObservablesMat <- sapply(1:d, function(j) {
  p <- pnorm(0.5 * x_true + rnorm(n, sd = 0.5))
  rbinom(n, 1, p)
})

## ----eval=TRUE----------------------------------------------------------------
library(lpmec)

# Run bootstrapped analysis
results <- lpmec(
  Y = Yobs,
  observables = as.data.frame(ObservablesMat),
  n_boot = 10,      # Reduced for demonstration
  n_partition = 5,  # Reduced for demonstration
  estimation_method = "em"
)

## ----eval=TRUE----------------------------------------------------------------
print(results)
summary(results)

## -----------------------------------------------------------------------------
plot(results)

## ----eval=TRUE----------------------------------------------------------------
# Bayesian MCMC estimation (requires Python environment setup)
if(FALSE){
mcmc_results <- lpmec(
  Y = Yobs,
  observables = as.data.frame(ObservablesMat),
  estimation_method = "mcmc",
  conda_env = "lpmec"  # Specify your conda environment
)
}

