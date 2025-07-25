% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpme_DoOneRun.R
\name{lpme_onerun}
\alias{lpme_onerun}
\title{lpme_onerun}
\usage{
lpme_onerun(
  Y,
  observables,
  observables_groupings = colnames(observables),
  make_observables_groupings = FALSE,
  estimation_method = "em",
  latent_estimation_fn = NULL,
  mcmc_control = list(backend = "pscl", n_samples_warmup = 500L, n_samples_mcmc = 1000L,
    batch_size = 512L, chain_method = "parallel", subsample_method = "full", n_thin_by =
    1L, n_chains = 2L),
  ordinal = FALSE,
  conda_env = "lpme",
  conda_env_required = TRUE
)
}
\arguments{
\item{Y}{A vector of observed outcome variables}

\item{observables}{A matrix of observable indicators used to estimate the latent variable}

\item{observables_groupings}{A vector specifying groupings for the observable indicators. Default is column names of observables.}

\item{make_observables_groupings}{Logical. If TRUE, creates dummy variables for each level of the observable indicators. Default is FALSE.}

\item{estimation_method}{Character specifying the estimation approach. Options include:
\itemize{
\item "em" (default): Uses expectation-maximization via \code{emIRT} package. Supports both binary (via \code{emIRT::binIRT}) and ordinal (via \code{emIRT::ordIRT}) indicators.
\item "averaging": Uses feature averaging.
\item "mcmc": Markov Chain Monte Carlo estimation using either \code{pscl::ideal} (R backend) or \code{numpyro} (Python backend)
\item "mcmc_joint": Joint Bayesian model that simultaneously estimates latent variables and outcome relationship using \code{numpyro}
\item "mcmc_overimputation": Two-stage MCMC approach with measurement error correction via over-imputation
\item "custom": In this case, latent estimation performed using \code{latent_estimation_fn}.
}}

\item{latent_estimation_fn}{Custom function for estimating latent trait from \code{observables} if \code{estimation_method="custom"} (optional). The function should accept a matrix of observables (rows are observations) and return a numeric vector of length equal to the number of observations.}

\item{mcmc_control}{A list indicating parameter specifications if MCMC used. 
\itemize{
  \item{\code{backend}}{
    Character string indicating the MCMC engine to use. Valid options are:
    \itemize{
      \item \code{"numpyro"} (default): Uses the Python \code{numpyro} package via \code{reticulate}.
      \item \code{"pscl"}: Uses the R-based \code{pscl::ideal} function.
    }
  }
  \item{\code{n_samples_warmup}}{
    Integer specifying the number of warm-up (a.k.a. burn-in) iterations
    before samples are collected. Default is \code{500}.
  }
  \item{\code{n_samples_mcmc}}{
    Integer specifying the number of post-warmup MCMC iterations to retain.
    Default is \code{1000}.
  }
  \item{\code{chain_method}}{
    Character string passed to \code{numpyro} specifying how to run multiple
    chains. Typical options include:
    \itemize{
      \item \code{"parallel"} (default): Runs chains in parallel.
      \item \code{"sequential"}: Runs chains sequentially.
      \item \code{"vectorized"}: Vectorized evaluation of multiple chains.
    }
  }
  \item{\code{n_thin_by}}{
    Integer indicating the thinning factor for MCMC samples (i.e., retaining
    every \code{n_thin_by}-th sample). Default is \code{1}.
  }
  \item{\code{n_chains}}{
    Integer specifying the number of parallel MCMC chains to run.
    Default is \code{2}.
  }
}}

\item{ordinal}{Logical indicating whether the observable indicators are ordinal (TRUE) or binary (FALSE).}

\item{conda_env}{A character string specifying the name of the conda environment to use 
via \code{reticulate}. Default is \code{"lpme"}.}

\item{conda_env_required}{A logical indicating whether the specified conda environment 
must be strictly used. If \code{TRUE}, an error is thrown if the environment is not found. 
Default is \code{TRUE}.}
}
\value{
A list containing various estimates and statistics:
\itemize{
  \item \code{ols_coef}: Coefficient from naive OLS regression
  \item \code{ols_se}: Standard error of naive OLS coefficient
  \item \code{ols_tstat}: T-statistic of naive OLS coefficient
  \item \code{corrected_ols_coef}: OLS coefficient corrected for measurement error
  \item \code{corrected_ols_se}: Standard error of corrected OLS coefficient (currently NA)
  \item \code{corrected_ols_tstat}: T-statistic of corrected OLS coefficient (currently NA)
  \item \code{corrected_ols_coef_alt}: Alternative corrected OLS coefficient (currently NA)
  \item \code{iv_coef}: Coefficient from instrumental variable regression
  \item \code{iv_se}: Standard error of IV regression coefficient
  \item \code{iv_tstat}: T-statistic of IV regression coefficient
  \item \code{corrected_iv_coef}: IV regression coefficient corrected for measurement error
  \item \code{corrected_iv_se}: Standard error of corrected IV coefficient
  \item \code{corrected_iv_tstat}: T-statistic of corrected IV coefficient
  \item \code{var_est_split}: Estimated variance of the measurement error
  \item \code{x_est1}: First set of latent variable estimates
  \item \code{x_est2}: Second set of latent variable estimates
}
}
\description{
Implements analysis for latent variable models with measurement error correction
}
\details{
This function implements a latent variable analysis with measurement error correction. 
It splits the observable indicators into two sets, estimates latent variables using each set, 
and then applies various correction methods including OLS correction and instrumental variable approaches.
}
\examples{
# Generate some example data
set.seed(123)
library( lpme )
Y <- rnorm(1000)
observables <- as.data.frame( matrix(sample(c(0,1), 1000*10, replace = TRUE), ncol = 10) )

# Run the analysis
results <- lpme_onerun(Y = Y, 
                       observables = observables)

# View the corrected estimates
print(results)

}
