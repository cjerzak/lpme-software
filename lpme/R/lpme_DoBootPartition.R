#' lpme
#'
#' Implements bootstrapped analysis for latent variable models with measurement error correction
#'
#' @param Y A vector of observed outcome variables
#' @param observables A matrix of observable indicators used to estimate the latent variable
#' @param orientation_signs (optional) A numeric vector of length equal to the number of columns in `observables`, containing 1 or -1 to indicate the desired orientation of each column. If provided, each column of `observables` will be oriented by this sign before analysis. Default is NULL (no orientation applied).
#' @param observables_groupings A vector specifying groupings for the observable indicators. Default is column names of observables.
#' @param make_observables_groupings Logical. If TRUE, creates dummy variables for each level of the observable indicators. Default is FALSE.
#' @param n_boot Integer. Number of bootstrap iterations. Default is 32.
#' @param n_partition Integer. Number of partitions for each bootstrap iteration. Default is 10.
#' @param boot_basis Vector of indices or grouping variable for stratified bootstrap. Default is 1:length(Y).
#' @param ordinal Logical indicating whether the observable indicators are ordinal (TRUE) or binary (FALSE).
#' @param return_intermediaries Logical. If TRUE, returns intermediate results. Default is TRUE.
#' @param estimation_method Character specifying the estimation approach. Options include:
#' \itemize{
#' \item "em" (default): Uses expectation-maximization via \code{emIRT} package. Supports both binary (via \code{emIRT::binIRT}) and ordinal (via \code{emIRT::ordIRT}) indicators.
#' \item "mcmc": Markov Chain Monte Carlo estimation using either \code{pscl::ideal} (R backend) or \code{numpyro} (Python backend)
#' \item "mcmc_joint": Full Bayesian model that simultaneously estimates latent variables and outcome relationship using \code{numpyro}
#' \item "mcmc_overimputation": Two-stage MCMC approach with measurement error correction via over-imputation
#' }
#' @param mcmc_control A list indicating parameter specifications if MCMC used. 
#' \itemize{
#'   \item{\code{backend}}{
#'     Character string indicating the MCMC engine to use. Valid options are:
#'     \itemize{
#'       \item \code{"numpyro"} (default): Uses the Python \code{numpyro} package via \code{reticulate}.
#'       \item \code{"pscl"}: Uses the R-based \code{pscl::ideal} function.
#'     }
#'   }
#'   \item{\code{n_samples_warmup}}{
#'     Integer specifying the number of warm-up (a.k.a. burn-in) iterations
#'     before samples are collected. Default is \code{500}.
#'   }
#'   \item{\code{n_samples_mcmc}}{
#'     Integer specifying the number of post-warmup MCMC iterations to retain.
#'     Default is \code{1000}.
#'   }
#'   \item{\code{chain_method}}{
#'     Character string passed to \code{numpyro} specifying how to run multiple
#'     chains. Typical options include:
#'     \itemize{
#'       \item \code{"parallel"} (default): Runs chains in parallel.
#'       \item \code{"sequential"}: Runs chains sequentially.
#'       \item \code{"vectorized"}: Vectorized evaluation of multiple chains.
#'     }
#'   }
#'   \item{\code{n_thin_by}}{
#'     Integer indicating the thinning factor for MCMC samples (i.e., retaining
#'     every \code{n_thin_by}-th sample). Default is \code{1}.
#'   }
#'   \item{\code{n_chains}}{
#'     Integer specifying the number of parallel MCMC chains to run.
#'     Default is \code{2}.
#'   }
#' }
#' @param conda_env A character string specifying the name of the conda environment to use 
#'   via \code{reticulate}. Default is \code{"lpme"}.
#' @param conda_env_required A logical indicating whether the specified conda environment 
#'   must be strictly used. If \code{TRUE}, an error is thrown if the environment is not found. 
#'   Default is \code{TRUE}.
#'
#' @return A list containing various estimates and statistics (in snake_case):
#' \itemize{
#'   \item \code{ols_coef}: Coefficient from naive OLS regression.
#'   \item \code{ols_se}: Standard error of naive OLS coefficient.
#'   \item \code{ols_tstat}: T-statistic of naive OLS coefficient.
#'   \item \code{iv_coef}: Coefficient from instrumental variable (IV) regression.
#'   \item \code{iv_se}: Standard error of IV regression coefficient.
#'   \item \code{iv_tstat}: T-statistic of IV regression coefficient.
#'   \item \code{corrected_iv_coef}: IV regression coefficient corrected for measurement error.
#'   \item \code{corrected_iv_se}: Standard error of the corrected IV coefficient (currently \code{NA}).
#'   \item \code{corrected_iv_tstat}: T-statistic of the corrected IV coefficient.
#'   \item \code{var_est}: Estimated variance of the measurement error (split-half variance).
#'   \item \code{corrected_ols_coef}: OLS coefficient corrected for measurement error.
#'   \item \code{corrected_ols_se}: Standard error of the corrected OLS coefficient (currently \code{NA}).
#'   \item \code{corrected_ols_tstat}: T-statistic of the corrected OLS coefficient (currently \code{NA}).
#'   \item \code{corrected_ols_coef_alt}: Alternative corrected OLS coefficient (if applicable).
#'   \item \code{corrected_ols_se_alt}: Standard error for the alternative corrected OLS coefficient (if applicable).
#'   \item \code{corrected_ols_tstat_alt}: T-statistic for the alternative corrected OLS coefficient (if applicable).
#'   \item \code{bayesian_ols_coef_outer_normed}: Posterior mean of the OLS coefficient under MCMC, 
#'     after normalizing by the overall sample standard deviation.
#'   \item \code{bayesian_ols_se_outer_normed}: Posterior standard error corresponding to \code{bayesian_ols_coef_outer_normed}.
#'   \item \code{bayesian_ols_tstat_outer_normed}: T-statistic for \code{bayesian_ols_coef_outer_normed}.
#'   \item \code{bayesian_ols_coef_inner_normed}: Posterior mean of the OLS coefficient under MCMC, 
#'     after normalizing each posterior draw individually.
#'   \item \code{bayesian_ols_se_inner_normed}: Posterior standard error corresponding to \code{bayesian_ols_coef_inner_normed}.
#'   \item \code{bayesian_ols_tstat_inner_normed}: T-statistic for \code{bayesian_ols_coef_inner_normed}.
#'   \item \code{m_stage_1_erv}: Extreme robustness value (ERV) for the first-stage regression 
#'     (\code{x_est2} on \code{x_est1}), if computed.
#'   \item \code{m_reduced_erv}: ERV for the reduced model (\code{Y} on \code{x_est1}), if computed.
#'   \item \code{x_est1}: First set of latent variable estimates.
#'   \item \code{x_est2}: Second set of latent variable estimates.
#' }
#'
#' @details 
#' This function implements a bootstrapped latent variable analysis with measurement error correction. 
#' It performs multiple bootstrap iterations, each with multiple partitions. For each partition, 
#' it calls the LatentOneRun function to estimate latent variables and apply various correction methods. 
#' The results are then aggregated across partitions and bootstrap iterations to produce final estimates 
#' and bootstrap standard errors.
#'
#' @examples
#' # Generate some example data
#' set.seed(123)
#' Y <- rnorm(1000)
#' observables <- as.data.frame( matrix(sample(c(0,1), 1000*10, replace = TRUE), ncol = 10) )
#' 
#' # Run the bootstrapped analysis
#' results <- lpme(Y = Y, 
#'                 observables = observables, 
#'                 n_boot = 10,    # small values for illustration only  
#'                 n_partition = 5 # small for size 
#'                 )
#' 
#' # View the corrected IV coefficient and its standard error
#' print(results)
#'
#' @export
#' @importFrom stats sd median 
#' @importFrom lpme lpme_onerun

lpme <- function(Y,
                 observables, 
                 observables_groupings = colnames(observables),
                 orientation_signs = NULL,
                 make_observables_groupings = FALSE,
                 n_boot = 32L, 
                 n_partition = 10L, 
                 boot_basis = 1:length(Y),
                 return_intermediaries = TRUE, 
                 ordinal = FALSE, 
                 estimation_method = "em",
                 mcmc_control = list(
                   backend = "numpyro",  # will override to use NumPyro-based MCMC
                   n_samples_warmup = 500L,
                   n_samples_mcmc   = 1000L,
                   batch_size = 512L, 
                   chain_method = "parallel", 
                   subsample_method = "full", 
                   anchor_parameter_id = NULL, 
                   n_thin_by = 1L, 
                   n_chains = 2L), 
                 conda_env = "lpme", 
                 conda_env_required = TRUE
                 ){ 
  
  # Orient the observables if orientation_signs are provided
  if (!is.null(orientation_signs)) {
    if (!is.numeric(orientation_signs) || length(orientation_signs) != ncol(observables)) {
      stop("orientation_signs must be a numeric vector of length equal to ncol(observables).")
    }
    if (!all(orientation_signs %in% c(1, -1))) {
      stop("orientation_signs must contain only 1 and -1.")
    }
    if(!all(unlist(observables) %in% c(0,1))){
      stop("Re-orientation in the non-binary case not yet implementated")
    }
    if(all(observables %in% c(0,1))){
      colnames_observables <- colnames(observables)
      observables <- sapply(1:ncol(observables), function(d_){
        observables[, d_] <- orientation_signs[d_] * observables[, d_] + 
                                    (1 - orientation_signs[d_]) / 2
      })
      colnames(observables) <- colnames_observables
    }
  }
  
  for(booti_ in seq_len(n_boot + 1L)){
    # if not the first iteration, sample bootstrap indices
    if(booti_ == 1L){
      boot_indices <- seq_along(Y)
    } else {
      # cluster/bootstrap over unique groups in boot_basis
      sampled_groups <- sample(
        unique(as.character(boot_basis)), 
        length(unique(boot_basis)), 
        replace = TRUE
      )
      # expand out to row indices
      boot_indices <- unlist(
        tapply(
          seq_along(boot_basis), 
          as.character(boot_basis), 
          c)[sampled_groups]
      )
    }
    
    for(parti_ in seq_len(n_partition)){
      message(sprintf("{booti_ %s of %s} -- {parti_ %s of %s}", booti_, n_boot+1, parti_, n_partition))
      
      # Run single analysis
      rungood <- F;runcounter <- 0; while(!rungood){ 
        runcounter <- runcounter + 1 
        LatentRunResults_ <- try(lpme_onerun(
          Y[boot_indices],
          observables[boot_indices,], 
          observables_groupings = observables_groupings,
          make_observables_groupings = make_observables_groupings, 
          estimation_method = estimation_method, 
          ordinal = ordinal, 
          mcmc_control = mcmc_control, 
          conda_env = conda_env,
          conda_env_required = conda_env_required
        ),T) 
        if(!"try-error" %in% class(LatentRunResults_)){ rungood <- TRUE }
        if(runcounter > 100){ stop("100 partition attempts failed... check data") }
      }
      
      
      # Tag each result with partition / bootstrap indices
      LatentRunResults_$PartitionIndex <- parti_
      LatentRunResults_$BootIndex      <- booti_
      
      # If first iteration, initialize the main results object
      if(booti_ == 1L && parti_ == 1L){
        LatentRunResults <- LatentRunResults_
      } else {
        # Otherwise, cbind new columns onto existing results
        for(name_ in names(LatentRunResults_)){
          LatentRunResults[[name_]] <- cbind( LatentRunResults[[name_]],  LatentRunResults_[[name_]] )
        }
      }
    }
  }
  
  # Summarizing function for across-partition/boot
  theSumFxn <- median
  # theSumFxn <- mean # (alternative if desired)
  
  # Now prepend "Intermediary_" to each piece of stored output
  names(LatentRunResults) <- paste0("Intermediary_", names(LatentRunResults))
  
  # Estimate cross-split variance of x_est1 - x_est2 for the first bootstrap
  # (Essentially the 'split-half' measure of measurement error)
  VarEst_split <- try(
    theSumFxn(
      apply(
        LatentRunResults$Intermediary_x_est1[, LatentRunResults$Intermediary_BootIndex == 1] -
          LatentRunResults$Intermediary_x_est2[, LatentRunResults$Intermediary_BootIndex == 1],
        1, sd
      )
    ),
    silent = TRUE
  )
  if(inherits(VarEst_split, "try-error")) { 
    VarEst_split <- NA 
  }
  
  # Standard error of the cross-split variance across bootstraps
  VarEst_split_se <- try(
    sd(
      sapply(2L:(n_boot + 1L), function(boot_){
        theSumFxn(
          apply(
            LatentRunResults$Intermediary_x_est1[, LatentRunResults$Intermediary_BootIndex == boot_] -
              LatentRunResults$Intermediary_x_est2[, LatentRunResults$Intermediary_BootIndex == boot_],
            1, sd
          )
        )
      })
    ),
    silent = TRUE
  )
  if(inherits(VarEst_split_se, "try-error")) { 
    VarEst_split_se <- NA 
  }
  
  # Helpers for summary stats
  qLow <- 0.025
  qUp  <- 0.975
  qf <- function(q, x) {
    if(length(x) < 2) return(NA_real_) 
    stats::quantile(x, prob = q, na.rm = TRUE)
  }
  
  # Pull out final estimates:
  # We apply tapply(..., theSumFxn) to the stored columns for each BootIndex,
  # then take the first element [1] of that vector (since we want the
  # median (or mean) across partitions for the "first" bootstrap, etc.).
  # Then we also produce an overall standard error across the bootstraps,
  # lower/upper bounds, etc.

  takeforse <- which(c(LatentRunResults$Intermediary_BootIndex)!=1)
  results <- list(
      # Naive OLS
      "ols_coef"   = (m1_ <- tapply(
        LatentRunResults$Intermediary_ols_coef, 
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "ols_se"     = (se1_ <- stats::sd(tapply(
        LatentRunResults$Intermediary_ols_coef[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      ))),
      "ols_lower"  = qf(qLow, tapply(
        LatentRunResults$Intermediary_ols_coef[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "ols_upper"  = qf(qUp, tapply(
        LatentRunResults$Intermediary_ols_coef[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "ols_tstat"  = (m1_ / se1_),
      
      # Corrected OLS
      "corrected_ols_coef_a" = (m1b_ <- tapply(
        LatentRunResults$Intermediary_corrected_ols_coef_a, 
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "corrected_ols_coef_b" = (m1b_ <- tapply(
        LatentRunResults$Intermediary_corrected_ols_coef_b, 
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "corrected_ols_coef" = (m1b_ <- tapply(
        LatentRunResults$Intermediary_corrected_ols_coef, 
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "corrected_ols_se"   = (se1b_ <- stats::sd(tapply(
        LatentRunResults$Intermediary_corrected_ols_coef[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      ))),
      "corrected_ols_lower" = qf(qLow, tapply(
        LatentRunResults$Intermediary_corrected_ols_coef[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "corrected_ols_upper" = qf(qUp, tapply(
        LatentRunResults$Intermediary_corrected_ols_coef[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "corrected_ols_tstat" = (m1b_ / se1b_),
      
      # Alternative corrected OLS
      "corrected_ols_coef_alt" = (m1b_ <- tapply(
        LatentRunResults$Intermediary_corrected_ols_coef_alt,
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "corrected_ols_se_alt"   = (se1b_ <- stats::sd(tapply(
        LatentRunResults$Intermediary_corrected_ols_coef_alt[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      ))),
      "corrected_ols_lower_alt" = qf(qLow, tapply(
        LatentRunResults$Intermediary_corrected_ols_coef_alt[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "corrected_ols_upper_alt" = qf(qUp, tapply(
        LatentRunResults$Intermediary_corrected_ols_coef_alt[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "corrected_ols_tstat_alt" = (m1b_ / se1b_),
      
      # IV regression
      "iv_coef_a" = (m2_ <- tapply(
        LatentRunResults$Intermediary_iv_coef_a, 
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "iv_coef_b" = (m2_ <- tapply(
        LatentRunResults$Intermediary_iv_coef_b, 
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "iv_coef" = (m2_ <- tapply(
        LatentRunResults$Intermediary_iv_coef, 
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "iv_se" = (se2_ <- stats::sd(tapply(
        LatentRunResults$Intermediary_iv_coef[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      ))),
      "iv_lower" = qf(qLow, tapply(
        LatentRunResults$Intermediary_iv_coef[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "iv_upper" = qf(qUp, tapply(
        LatentRunResults$Intermediary_iv_coef[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "iv_tstat" = (m2_ / se2_),
      
      # Corrected IV
      "corrected_iv_coef_a" = (m2_ <- tapply(
        LatentRunResults$Intermediary_corrected_iv_coef_a, 
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "corrected_iv_coef_b" = (m2_ <- tapply(
        LatentRunResults$Intermediary_corrected_iv_coef_b, 
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "corrected_iv_coef" = (m4_ <- tapply(
        LatentRunResults$Intermediary_corrected_iv_coef, 
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "corrected_iv_se" = (se4_ <- stats::sd(tapply(
        LatentRunResults$Intermediary_corrected_iv_coef[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      ))),
      "corrected_iv_lower" = qf(qLow, tapply(
        LatentRunResults$Intermediary_corrected_iv_coef[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "corrected_iv_upper" = qf(qUp, tapply(
        LatentRunResults$Intermediary_corrected_iv_coef[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "corrected_iv_tstat" = (m4_ / se4_),
      
      # Bayesian OLS (outer-normed)
      "bayesian_ols_coef_outer_normed" = (m4_ <- tapply(
        LatentRunResults$Intermediary_bayesian_ols_coef_outer_normed,
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "bayesian_ols_se_outer_normed_parametric" = (se4_b <- tapply(
        LatentRunResults$Intermediary_bayesian_ols_se_outer_normed,
        LatentRunResults$Intermediary_BootIndex,
        function(x) { 1/length(x) * sqrt(sum(x^2)) }
      )[1]),
      "bayesian_ols_se_outer_normed" = (se4_ <- stats::sd(tapply(
        LatentRunResults$Intermediary_bayesian_ols_coef_outer_normed[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      ))),
      "bayesian_ols_lower_outer_normed" = qf(qLow, tapply(
        LatentRunResults$Intermediary_bayesian_ols_coef_outer_normed[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "bayesian_ols_upper_outer_normed" = qf(qUp, tapply(
        LatentRunResults$Intermediary_bayesian_ols_coef_outer_normed[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "bayesian_ols_tstat_outer_normed" = (m4_ / se4_),
      
      # Bayesian OLS (inner-normed)
      "bayesian_ols_coef_inner_normed" = (m4_ <- tapply(
        LatentRunResults$Intermediary_bayesian_ols_coef_inner_normed, 
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "bayesian_ols_se_inner_normed_parametric" = (se4_b <- tapply(
        LatentRunResults$Intermediary_bayesian_ols_se_inner_normed[takeforse],
        LatentRunResults$Intermediary_BootIndex[takeforse],
        function(x) { 1/length(x) * sqrt(sum(x^2)) }
      )[1]),
      "bayesian_ols_se_inner_normed" = (se4_ <- stats::sd(tapply(
        LatentRunResults$Intermediary_bayesian_ols_coef_inner_normed[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      ))),
      "bayesian_ols_lower_inner_normed" = qf(qLow, tapply(
        LatentRunResults$Intermediary_bayesian_ols_coef_inner_normed[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "bayesian_ols_upper_inner_normed" = qf(qUp, tapply(
        LatentRunResults$Intermediary_bayesian_ols_coef_inner_normed[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "bayesian_ols_tstat_inner_normed" = (m4_ / se4_),
      
      # Robustness-value measures 
      "m_stage_1_erv" = (m2_ <- tapply(
        LatentRunResults$Intermediary_m_stage_1_erv, 
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "m_stage_1_erv_se" = (se2_ <- stats::sd(tapply(
        LatentRunResults$Intermediary_m_stage_1_erv[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      ))),
      "m_stage_1_erv_lower" = qf(qLow, tapply(
        LatentRunResults$Intermediary_m_stage_1_erv[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "m_stage_1_erv_upper" = qf(qUp, tapply(
        LatentRunResults$Intermediary_m_stage_1_erv[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "m_stage_1_erv_tstat" = (m2_ / se2_),
      
      "m_reduced_erv" = (m2_ <- tapply(
        LatentRunResults$Intermediary_m_reduced_erv, 
        LatentRunResults$Intermediary_BootIndex, theSumFxn
      )[1]),
      "m_reduced_erv_se" = (se2_ <- stats::sd(tapply(
        LatentRunResults$Intermediary_m_reduced_erv[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      ))),
      "m_reduced_erv_lower" = qf(qLow, tapply(
        LatentRunResults$Intermediary_m_reduced_erv[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "m_reduced_erv_upper" = qf(qUp, tapply(
        LatentRunResults$Intermediary_m_reduced_erv[takeforse], 
        LatentRunResults$Intermediary_BootIndex[takeforse], theSumFxn
      )),
      "m_reduced_erv_tstat" = (m2_ / se2_),
      
      # Final single-run estimates for x
      "x_est1" = LatentRunResults$Intermediary_x_est1[, 1],
      "x_est2" = LatentRunResults$Intermediary_x_est2[, 1],
      
      # Additional summary of the cross-split variance
      "var_est_split"    = VarEst_split,
      "var_est_split_se" = VarEst_split_se
    )
  class(results) <- "lpme"
  return(results)
  
}
