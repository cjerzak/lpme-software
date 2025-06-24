#' lpme_onerun
#'
#' Implements analysis for latent variable models with measurement error correction
#'
#' @param Y A vector of observed outcome variables
#' @param observables A matrix of observable indicators used to estimate the latent variable
#' @param observables_groupings A vector specifying groupings for the observable indicators. Default is column names of observables.
#' @param make_observables_groupings Logical. If TRUE, creates dummy variables for each level of the observable indicators. Default is FALSE.
#' @param estimation_method Character specifying the estimation approach. Options include:
#' \itemize{
#' \item "em" (default): Uses expectation-maximization via \code{emIRT} package. Supports both binary (via \code{emIRT::binIRT}) and ordinal (via \code{emIRT::ordIRT}) indicators.
#' \item "averaging": Uses feature averaging.
#' \item "mcmc": Markov Chain Monte Carlo estimation using either \code{pscl::ideal} (R backend) or \code{numpyro} (Python backend)
#' \item "mcmc_joint": Joint Bayesian model that simultaneously estimates latent variables and outcome relationship using \code{numpyro}
#' \item "mcmc_overimputation": Two-stage MCMC approach with measurement error correction via over-imputation
#' \item "custom": In this case, latent estimation performed using \code{latent_estimation_fn}.
#' }
#' @param latent_estimation_fn Custom function for estimating latent trait from \code{observables} if \code{estimation_method="custom"} (optional). The function should accept a matrix of observables (rows are observations) and return a numeric vector of length equal to the number of observations.
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
#' @param ordinal Logical indicating whether the observable indicators are ordinal (TRUE) or binary (FALSE).
#'
#' @return A list containing various estimates and statistics:
#' \itemize{
#'   \item \code{ols_coef}: Coefficient from naive OLS regression
#'   \item \code{ols_se}: Standard error of naive OLS coefficient
#'   \item \code{ols_tstat}: T-statistic of naive OLS coefficient
#'   \item \code{corrected_ols_coef}: OLS coefficient corrected for measurement error
#'   \item \code{corrected_ols_se}: Standard error of corrected OLS coefficient (currently NA)
#'   \item \code{corrected_ols_tstat}: T-statistic of corrected OLS coefficient (currently NA)
#'   \item \code{corrected_ols_coef_alt}: Alternative corrected OLS coefficient (currently NA)
#'   \item \code{iv_coef}: Coefficient from instrumental variable regression
#'   \item \code{iv_se}: Standard error of IV regression coefficient
#'   \item \code{iv_tstat}: T-statistic of IV regression coefficient
#'   \item \code{corrected_iv_coef}: IV regression coefficient corrected for measurement error
#'   \item \code{corrected_iv_se}: Standard error of corrected IV coefficient
#'   \item \code{corrected_iv_tstat}: T-statistic of corrected IV coefficient
#'   \item \code{var_est_split}: Estimated variance of the measurement error
#'   \item \code{x_est1}: First set of latent variable estimates
#'   \item \code{x_est2}: Second set of latent variable estimates
#' }
#'
#' @details 
#' This function implements a latent variable analysis with measurement error correction. 
#' It splits the observable indicators into two sets, estimates latent variables using each set, 
#' and then applies various correction methods including OLS correction and instrumental variable approaches.
#'
#'
#' @examples
#' # Generate some example data
#' set.seed(123)
#' library( lpme )
#' Y <- rnorm(1000)
#' observables <- as.data.frame( matrix(sample(c(0,1), 1000*10, replace = TRUE), ncol = 10) )
#' 
#' # Run the analysis
#' results <- lpme_onerun(Y = Y, 
#'                        observables = observables)
#' 
#' # View the corrected estimates
#' print(results)
#'
#' @export
#' @importFrom stats lm cor var rnorm
#' @importFrom AER ivreg
#' @importFrom pscl rollcall
#' @importFrom sandwich vcovHC
#' @importFrom mvtnorm rmvnorm
#' @importFrom Amelia amelia
#' @importFrom gtools quantcut
#' @importFrom emIRT binIRT ordIRT
#' @importFrom sensemakr extreme_robustness_value

lpme_onerun <- function( Y,
                         observables, 
                         observables_groupings = colnames(observables),
                         make_observables_groupings = FALSE, 
                         estimation_method = "em", 
                         mcmc_control = list(
                           backend = "numpyro",  # will override to use NumPyro-based MCMC
                           n_samples_warmup = 500L,
                           n_samples_mcmc   = 1000L,
                           batch_size = 512L, 
                           chain_method = "parallel", 
                           subsample_method = "full", 
                           n_thin_by = 1L, 
                           n_chains = 2L), 
                         ordinal = FALSE, 
                         conda_env = "lpme", 
                         conda_env_required = TRUE){
  # coerce to data.frame
  observables <- as.data.frame( observables )
  
  t0 <- Sys.time()
  INIT_SCALER <- 1/10
  Bayesian_OLSSE_InnerNormed <- Bayesian_OLSCoef_InnerNormed <- NA; 
  Bayesian_OLSSE_OuterNormed <- Bayesian_OLSCoef_OuterNormed <- NA;
  FullBayesianSlope_mean <- FullBayesianSlope_std <- NA
  items.split1_names <- sample(unique(observables_groupings), 
                               size = floor(length(unique(observables_groupings))/2), replace=FALSE)
  items.split2_names <- unique(observables_groupings)[! (observables_groupings %in% items.split1_names)]
  for(split_ in c("", "1", "2")){
    if(split_ == ""){ items.split_ <- 1:length(observables_groupings) }
    if(split_ == "1"){ items.split_ <- (1:length(observables_groupings))[observables_groupings %in% items.split1_names] }
    if(split_ == "2"){ items.split_ <- (1:length(observables_groupings))[observables_groupings %in% items.split2_names] }
    
    # estimating ideal points
    if(make_observables_groupings == FALSE){
      observables_ <- observables[,items.split_]
    }
    if(make_observables_groupings == TRUE){
      observables_ <- do.call(cbind,unlist(apply(observables[,items.split_],2,function(zer){
        list( model.matrix(~0+as.factor(zer))[,-1] )}),recursive = FALSE))
    }
    
    if(estimation_method == "pca"){
        # For PCA, we typically want numeric inputs only
        x_init <- scale(apply( observables_, 1, function(x){ mean(f2n(x), na.rm=TRUE)})) * INIT_SCALER
        observables__ <- apply(as.matrix(observables_), 2, f2n)
        
        # do not zero impute
        #observables__[is.na(observables__)] <- 0
        
        # Run PCA on centered + scaled columns:
        pca_out <- prcomp(observables__, center = TRUE, scale. = TRUE)
        
        # Extract the first principal component
        x.est_ <- pca_out$x[, 1]
        
        # Scale so that it has mean 0, sd 1
        x.est_ <- scale(x.est_)
        
        # (Optional) Flip sign to match an anchor — 
        # compare with whichever initial reference you want, e.g. x_init
        if (x.est_[which.max(x_init)] < 0) {
          x.est_ <- -1 * x.est_
        }
    }
    
    if (estimation_method == "averaging") {
      # Final estimate is just this averaged measure (split-by-split)
      x.est_ <- scale(apply(observables_, 1, function(x) mean(f2n(x), na.rm = TRUE)))
    }
    
    if (estimation_method == "custom") {
      # Final estimate is just this averaged measure (split-by-split)
      x.est_ <- scale( latent_estimation_fn(observables_) )
    }
    
    # first, do emIRT
    if(estimation_method == "em"){
      x_init <- scale(apply( observables_, 1, function(x){ mean(f2n(x), na.rm=TRUE)})) * INIT_SCALER

      if(ordinal){
        observables__ <- apply(as.matrix(observables_),2,f2n)
        if( !all(c(observables__) %in% 1:3) ){
          observables__ <- apply(observables__, 2, function(x){ 
            r <- rank(f2n(x), ties.method = "average")
            group <- ceiling(r / (length(x) / 3))
            return(group)
          })
        }
        if( !all(apply(observables__,2,function(x)length(unique(x))) == 3) ){
          warning("Some observables mapped to binary, not ordinal values -- dropping those.")
          observables__ <- observables__[,which(apply(observables__,2,function(x){length(unique(x))})==3)]
          # apply(observables__,2,table)
          # colSums(apply(observables__,2,table))
          # sum(observables__)
          # summary(c(cor(observables__)[lower.tri(cor(observables__))]))
        }
        observables__[is.na(observables__)] <- 0 # map missing to zero
        
        ## Generate starts and priors for ordered model
        JJ <- ncol(observables__)
        NN <- nrow(observables__)
        
        # starts
        starts <- vector(mode = "list")
        starts$DD <- matrix(rep(0.5,JJ), ncol=1)
        starts$tau <- matrix(rep(-0.5,JJ), ncol=1)
        starts$beta <- matrix(runif(JJ,-1,1), ncol=1) 
        starts$x <- as.matrix(c(x_init))
        
        priors <- list("x" = list(mu = matrix(0,1,1), sigma = matrix(1,1,1) ), 
                      "beta" = list(mu = matrix(0,2,1), sigma = matrix(diag(25,2),2,2)))
        
        capture.output(
          out_emIRT <- try(emIRT::ordIRT(.rc = observables__,
                                  .starts = starts, 
                                  .priors = priors, 
                                  .control = list(
                                    threads = 1,
                                    verbose = FALSE,
                                    thresh = 1e-6,
                                    maxit=3000,
                                    checkfreq=50)
                        ),TRUE)
        )
        if("try-error" %in% class(out_emIRT)){ 
          stop(out_emIRT)
        } 
        
        # Scale the final ideal points and store
        x.est_ <- scale(out_emIRT$means$x)
        if(x.est_[which.max(x_init)] < 0){ x.est_ <- -1*x.est_ }  # anchor - no .anchor_subject arg 
        
        # FORCING FORCING SANITY SANITY
        #x.est_ <- scale(x_init)
        
        x.est_EM <- x.est_
      }
    
      if( !ordinal ){ 
        # informative initialization 
        rc_ <- emIRT::convertRC( pscl::rollcall(observables_) )
        capture.output( out_emIRT <- emIRT::binIRT(.rc = rc_, 
                            .starts = list("alpha" = matrix(rnorm(ncol(observables_), sd = .1)),
                                           "beta" = matrix(rnorm(ncol(observables_), sd = 1)),
                                           "x" = matrix(x_init)), 
                            .priors = emIRT::makePriors(.N= rc_$n, .J = rc_$m, .D = 1), 
                            .control= list(threads=1, verbose=FALSE, thresh=1e-6,verbose=FALSE),
                            .anchor_subject = which.max(x_init)
                            ) ) 
        x.est_EM <- x.est_ <- scale(out_emIRT$means$x); 
      }
    }
    
    if( grepl(estimation_method, pattern = "mcmc") & mcmc_control$backend == "pscl" ){
      if(estimation_method == "mcmc_joint"){stop('Must use numpyro with option: estimation_method="mcmc_joint"')}
      if(estimation_method == "mcmc" ){
        t0_ <- Sys.time()
        startval_ <- rowMeans( observables_ )
        maxiter_ <- mcmc_control$n_samples_mcmc + mcmc_control$n_samples_warmup
        burnin_ <- mcmc_control$n_samples_warmup
        thin_ <- mcmc_control$n_thin_by
        pscl_input_ <- pscl::rollcall(observables_)
        capture.output( 
          pscl_ideal <- pscl::ideal( pscl_input_, 
                            normalize = TRUE, 
                            store.item = TRUE, 
                            startvals = list("x" = startval_ ),
                            maxiter = maxiter_, 
                            burnin = burnin_,
                            thin = thin_)
        )
        message(sprintf("\n MCMC Runtime: %.3f min",  tdiff_ <- as.numeric(difftime(Sys.time(),  t0_, units = "secs"))/60))
        # mean(pscl_ideal$xbar); sd(pscl_ideal$xbar) # confirm 0 and 1 
        x.est_MCMC <- x.est_ <- pscl_ideal$xbar; s_past <- 1 # summary(lm(Y~x.est_))
        if( split_ == "" ){
          RescaledAbilities_OuterNormed  <- t(pscl_ideal$x[,,1])
          #Bayesian_OLSCoef_OuterNormed <- apply(RescaledAbilities_OuterNormed, 2, function(x_){ coef(lm(Y~x_))[2]}) # simple MOC 
          Bayesian_OLSCoef_OuterNormed <- apply(RescaledAbilities_OuterNormed, 2, function(x_){ 
            VCovHat <- sandwich::vcovHC( myModel <- lm(Y~x_), type = "HC3" ) 
            coef_ <- mvtnorm::rmvnorm(n = 1, mean = coef(myModel), sigma = VCovHat)[1,2]
          } ) # complicated MOC 
          
          # summary( apply(RescaledAbilities_outer, 2, sd) ) 
          # sd(rowMeans(RescaledAbilities_outer));mean(rowMeans(RescaledAbilities_outer)) # confirm sanity value of 1, 0
          # hist( Bayesian_OLSCoef_OuterNormed ); summary( Bayesian_OLSCoef_OuterNormed)
          Bayesian_OLSSE_OuterNormed <- sd( Bayesian_OLSCoef_OuterNormed ) 
          Bayesian_OLSCoef_OuterNormed <- mean( Bayesian_OLSCoef_OuterNormed )
          
          # InnerNormed
          RescaledAbilities_InnerNormed  <- ( apply(t(pscl_ideal$x[,,1]), 2, function(x_){scale(x_)})  ) 
          Bayesian_OLSCoef_InnerNormed <- apply(RescaledAbilities_InnerNormed, 2, function(x_){ 
            VCovHat <- sandwich::vcovHC( myModel <- lm(Y~x_), type = "HC3" ) 
            coef_ <- mvtnorm::rmvnorm(n = 1, mean = coef(myModel), sigma = VCovHat)[1,2]
          } ) # complicated MOC 
          # apply(RescaledAbilities, 2, sd)
          # sd(rowMeans(RescaledAbilities));mean(rowMeans(RescaledAbilities)) # confirm sanity value of 1, 0
          # hist( Bayesian_OLSCoef_InnerNormed ); summary( Bayesian_OLSCoef_InnerNormed )
          Bayesian_OLSSE_InnerNormed <- sd( Bayesian_OLSCoef_InnerNormed ) 
          Bayesian_OLSCoef_InnerNormed <- mean( Bayesian_OLSCoef_InnerNormed )
        }
      }
      if(estimation_method == "mcmc_overimputation" & split_ == ""){
        t0_ <- Sys.time()
        startval_ <- rowMeans( observables_ )
        maxiter_ <- mcmc_control$n_samples_mcmc + mcmc_control$n_samples_warmup
        burnin_ <- mcmc_control$n_samples_warmup
        thin_ <- mcmc_control$n_thin_by
        pscl_input_ <- pscl::rollcall(observables_)
        capture.output( 
          pscl_ideal <- pscl::ideal( pscl_input_, 
                                     normalize = TRUE, 
                                     store.item = TRUE, 
                                     startvals = list("x" = startval_ ),
                                     maxiter = maxiter_, 
                                     burnin = burnin_,
                                     thin = thin_)
        )
        # mean(pscl_ideal$xbar); sd(pscl_ideal$xbar) # confirm 0 and 1 
        x.est_MCMC <- x.est_ <- pscl_ideal$xbar; s_past <- 1 
        if( split_ == "" ){
            for(outType_ in c("Outer","Inner")){ 
              # file:///Users/cjerzak/Dropbox/LatentMeasures/literature/CAUGHEY-ps8-solution.html
              if(outType_ == "Outer"){
                RescaledAbilities  <- t(pscl_ideal$x[,,1])
              }
              if(outType_ == "Inner"){
                RescaledAbilities  <- t(pscl_ideal$x[,,1])
                RescaledAbilities  <- ( apply(RescaledAbilities, 2, function(x_){scale(x_)})  ) 
              }
              Xobs_mean <- apply(RescaledAbilities, 1, function(x_){ mean(x_) } ) 
              Xobs_SE <- apply(RescaledAbilities, 1, function(x_){ sd(x_) } ) 
              
              dat_ <- cbind(Y, x.est_MCMC)
              
              # Specify (over-)imputation model
              # priors: #a numeric matrix with four columns 
              # (row, column, mean, standard deviation) 
              # indicating noisy estimates of the values to be imputed.
              outcome_priors <- cbind(
                1:nrow(dat_), 1,              
                Y, # mean 
                1 # sd 
              )
              policy_priors <- cbind( 
                1:nrow(dat_), 2,              
                Xobs_mean,
                Xobs_SE
              )
              #priors <- rbind(policy_priors, outcome_priors)
              priors <- policy_priors
              
              #a numeric matrix where each row indicates a row and column of x to be overimputed.
              # overimp <- rbind(cbind(1:nrow(dat_), 1), cbind(1:nrow(dat_), 2))
              overimp <- cbind(1:nrow(dat_), 2)
              
              # perform overimputation 
              overimputed_data <- Amelia::amelia( 
                x = dat_, 
                m = (nOverImpute <- 5), # default 
                p2s = 0,
                priors = priors, 
                overimp = overimp, 
                parallel = "no")
              
              # Perform multiple overimputation
              overimputed_Y <- do.call(cbind,lapply(overimputed_data$imputations,function(l_){l_[,1]}))
              overimputed_x.est_MCMC <- do.call(cbind,lapply(overimputed_data$imputations,function(l_){l_[,2]}))
              if(outType_ == "Outer"){
                overimputed_x.est_MCMC  <- (overimputed_x.est_MCMC-mean(rowMeans(overimputed_x.est_MCMC)))/sd(rowMeans(overimputed_x.est_MCMC))
                #sd(rowMeans(overimputed_x.est_MCMC))
              }
              if(outType_ == "Inner"){
                overimputed_x.est_MCMC  <- ( apply(overimputed_x.est_MCMC, 2, function(x_){scale(x_)})  ) 
                #apply(overimputed_x.est_MCMC,2,sd)
              }
              # cor(cbind(x.est_MCMC, overimputed_x.est_MCMC))

              # Analyze overimputed datasets
              overimputed_coefs <- unlist(unlist( sapply(1:nOverImpute, function(s_){
                coef(lm(overimputed_Y[,s_] ~ overimputed_x.est_MCMC[,s_]))[2]
              })))
              if(outType_ == "Outer"){ 
                Bayesian_OLSSE_OuterNormed <- sd( overimputed_coefs )
                Bayesian_OLSCoef_OuterNormed <- mean( overimputed_coefs )
              }
              if(outType_ == "Inner"){ 
                Bayesian_OLSSE_InnerNormed <- sd( overimputed_coefs )
                Bayesian_OLSCoef_InnerNormed <- mean( overimputed_coefs )
              }
            }
        }
        message(sprintf("\n Overimputation Runtime: %.3f min",  tdiff_ <- as.numeric(difftime(Sys.time(),  t0_, units = "secs"))/60))
      }
      if(estimation_method == "mcmc_joint" & split_ == ""){
        ## ?? 
      }
    }
    if( grepl(estimation_method, pattern = "mcmc") & mcmc_control$backend == "numpyro" ){
        if(!"jax" %in% ls(envir = lpme_env)){ initialize_jax(conda_env, conda_env_required) }
        
        # Construct for annotating conditionally independent variables.
        # Within a plate context manager, sample sites will be automatically broadcasted to the size of the plate. 
        # Additionally, a scale factor might be applied by certain inference algorithms
        # ifsubsample_size is specified.
        
        # Set up MCMC
        lpme_env$numpyro$set_host_device_count( mcmc_control$n_chains )
        N <- ai(nrow(observables_))
        K <- ai(ncol(observables_))
        
        # Define the two-parameter IRT model using Matt's trick + subsampling
        # Note: IRTModel_batch is depreciated 
        IRTModel_batch <- function(X, Y) {
            # Number of observations (rows) and items (columns)
            N <- X$shape[[1]]
            K <- X$shape[[2]]
            
            # Global hyperpriors for ability (non-centered)
            mu_ability <- lpme_env$numpyro$sample("mu_ability",
                                                  lpme_env$dist$Normal(0, 1))
            sigma_ability <- lpme_env$numpyro$sample("sigma_ability",
                                                     lpme_env$dist$HalfNormal(1))
            with(lpme_env$numpyro$plate("rows", N), {
              eps_ability <- lpme_env$numpyro$sample("eps_ability",
                                                     lpme_env$dist$Normal(0, 1))
              ability <- lpme_env$numpyro$deterministic(
                "ability", (mu_ability + sigma_ability * eps_ability)[,NULL] )
            })
            
            # Hyperpriors for item parameters
            mu_difficulty <- lpme_env$numpyro$sample("mu_difficulty",
                                                     lpme_env$dist$Normal(0, 3))
            sigma_difficulty <- lpme_env$numpyro$sample("sigma_difficulty",
                                                        lpme_env$dist$HalfNormal(3))
            mu_log_discrimination <- lpme_env$numpyro$sample("mu_log_discrimination",
                                                             lpme_env$dist$Normal(0.5, 1))
            sigma_log_discrimination <- lpme_env$numpyro$sample("sigma_log_discrimination",
                                                                lpme_env$dist$HalfNormal(0.5))
            with(lpme_env$numpyro$plate("columns", K), {
              eps_difficulty <- lpme_env$numpyro$sample("eps_difficulty",
                                                        lpme_env$dist$Normal(0, 3))
              difficulty <- lpme_env$numpyro$deterministic(
                "difficulty", (mu_difficulty + sigma_difficulty * eps_difficulty)[NULL,] )
              
              eps_log_discrimination <- lpme_env$numpyro$sample("eps_log_discrimination",
                                                      lpme_env$dist$Normal(0, 1))
              log_discrimination <- mu_log_discrimination + sigma_log_discrimination * eps_log_discrimination
              discrimination <- lpme_env$numpyro$deterministic(
                "discrimination", lpme_env$jax$nn$softplus(log_discrimination)[NULL,] )
            })
            
            # Define a local likelihood function for the *subset* of rows
            local_lik_fn <- function(){
              with(lpme_env$numpyro$plate(
                "rows_subsample",
                size = N,
                subsample_size = ai(mcmc_control$batch_size),
                dim = -2L
              ) %as% "idx", {
                # Subset 
                ability_sub <- lpme_env$jnp$take(ability, idx, axis = 0L)
                
                # Observed subset of X
                X_sub <- lpme_env$jnp$take(X, indices = idx, axis = 0L)
                
                # Construct logit for subset
                logits_sub <- (ability_sub - difficulty) * discrimination
                
                # Add columns plate around the sample statement
                with(lpme_env$numpyro$plate("columns", K, dim = -1L), {
                  lpme_env$numpyro$sample(
                    "Xlik_sub",
                    lpme_env$dist$Bernoulli(logits = logits_sub),
                    obs = X_sub )
                })
                
                # If doing a full regression on Y:
                if (estimation_method == "mcmc_joint") {
                  Y_intercept <- lpme_env$numpyro$sample("YModel_intercept",
                                                         lpme_env$dist$Normal(0, 1))
                  Y_slope <- lpme_env$numpyro$sample("YModel_slope",
                                                     lpme_env$dist$Normal(0, 1))
                  Y_sigma <- lpme_env$numpyro$sample("YModel_sigma",
                                                     lpme_env$dist$HalfNormal(1))
                  
                  Y_sub <- lpme_env$jnp$take(Y, idx)
                  Y_mu_sub <- Y_intercept + Y_slope * ability_sub
                  lpme_env$numpyro$sample(
                    "Ylik_sub",
                    lpme_env$dist$Normal(Y_mu_sub, Y_sigma),
                    obs = Y_sub
                  )
                }
              })
            }
            
            # Scale the local likelihood by (N / batch_size) for unbiased gradient
            scaled_lik <- lpme_env$numpyro$handlers$scale(
              scale = as.numeric(N / mcmc_control$batch_size))(local_lik_fn)
            
            # Execute the scaled likelihood
            #scaled_lik() # documentation doesn't seem to scale?
            local_lik_fn()
          }
        
        # Define the two-parameter IRT model using Matt's trick 
        IRTModel_full <- (function(X, # binary indicators 
                                   Y  # outcome 
        ){
            
            # Global hyperpriors for ability
            mu_ability <- lpme_env$numpyro$sample("mu_ability",
                                                  lpme_env$dist$Normal(0, 1))
            sigma_ability <- lpme_env$numpyro$sample("sigma_ability",
                                                     lpme_env$dist$HalfNormal(0.5))
            
            # Non-centered parameterization for ability
            with(lpme_env$numpyro$plate("rows", X$shape[[1]], dim = -2L), { # N
              eps_ability <- lpme_env$numpyro$sample("eps_ability",
                                                     lpme_env$dist$Normal(0, 1))
              ability <- lpme_env$numpyro$deterministic("ability", 
                                 (mu_ability + sigma_ability * eps_ability) )
            })
            
            # If the user gave us an anchor_parameter_id, we convert it to 0-based
            # for Python indexing. We will forcibly ensure difficulty[anchor_id] is > 0.
            anchor_id <- NULL
            if (!is.null(mcmc_control$anchor_parameter_id)) {
              anchor_id <- as.integer(mcmc_control$anchor_parameter_id - 1L)  # R is 1-based, NumPyro is 0-based
            }
            
            # For difficulty, you might keep it simple or also do a non-centered param:
            mu_difficulty <- lpme_env$numpyro$sample("mu_difficulty",
                                                     lpme_env$dist$Normal(0, 2))
            sigma_difficulty <- lpme_env$numpyro$sample("sigma_difficulty",
                                                        lpme_env$dist$HalfNormal(1))
            mu_log_discrimination <- lpme_env$numpyro$sample("mu_log_discrimination",
                                                             lpme_env$dist$Normal(0.5, 1))
            sigma_log_discrimination <- lpme_env$numpyro$sample("sigma_log_discrimination",
                                                                lpme_env$dist$HalfNormal(0.5))
            with(lpme_env$numpyro$plate("columns", X$shape[[2]], dim = -1L), {  # D
              difficulty_raw <- lpme_env$numpyro$sample("eps_difficulty",
                                                        lpme_env$dist$Normal(0, 1))
              
              # If anchor_id was specified, force difficulty for that item to be positive
              difficulty_adjusted <- difficulty_raw
              if( !is.null(anchor_id) ){
                # Use at update to ensure the anchor difficulty is strictly positive
                difficulty_adjusted <- difficulty_raw$at[anchor_id]$set(
                  lpme_env$jax$nn$softplus(difficulty_raw[anchor_id])
                )
              }
              difficulty <- lpme_env$numpyro$deterministic("difficulty", 
                                         (difficulty_adjusted)[NULL,]) # if NOT using centered parameterization
                                        #(mu_difficulty + sigma_difficulty * difficulty_adjusted)[NULL,]) # if using centered parameterization

              # discrimination <- lpme_env$numpyro$sample("discrimination", lpme_env$dist$HalfNormal(2))
              eps_log_discrimination <- lpme_env$numpyro$sample("eps_log_discrimination", 
                                                                lpme_env$dist$Normal(0.5, 2))
              discrimination <- lpme_env$numpyro$deterministic("discrimination", 
                                        lpme_env$jax$nn$softplus(eps_log_discrimination)[NULL,]) # if NOT using centered parameterization
                                        #lpme_env$jax$nn$softplus(mu_log_discrimination + sigma_log_discrimination * eps_log_discrimination)[NULL,]) # if using centered parameterization
            })
            
            # Construct logits using the new 'ability' & 'difficulty'
            # note: discrimination constrained to be positive 
            Xprobs <- ( ability - difficulty ) * discrimination
            
            # Likelihood for X using a Bernoulli logit, wrapped in plates
            #lpme_env$numpyro$sample("Xlik", lpme_env$dist$Bernoulli(logits=Xprobs), obs=X)
            with(lpme_env$numpyro$plate("rows", X$shape[[1]], dim = -2L), {
              with(lpme_env$numpyro$plate("columns", X$shape[[2]], dim = -1L), {
                lpme_env$numpyro$sample("Xlik", lpme_env$dist$Bernoulli(logits = Xprobs), obs = X)
              }) })
            
            # If you are using the outcome portion in "mcmc_joint" mode:
            if(estimation_method == "mcmc_joint"){
              Y_intercept <- lpme_env$numpyro$sample("YModel_intercept",
                                                     lpme_env$dist$Normal(0, 1))
              Y_slope <- lpme_env$numpyro$sample("YModel_slope",
                                                 lpme_env$dist$Normal(0, 1))
              Y_sigma <- lpme_env$numpyro$sample("YModel_sigma",
                                                 lpme_env$dist$HalfNormal(1))
              
              Y_mu <- Y_intercept + Y_slope * ability
              lpme_env$numpyro$sample("Ylik",
                                      lpme_env$dist$Normal(Y_mu, Y_sigma),
                                      obs=Y)
            }
          })
        
      # setup & run a MCMC run
      if(mcmc_control$subsample_method == "batch"){ 
        message("Enlisting HMCECS kernels...")
        # https://num.pyro.ai/en/stable/mcmc.html#numpyro.infer.hmc_gibbs.HMCECS
        kernel <- lpme_env$numpyro$infer$HMCECS(lpme_env$numpyro$infer$NUTS(IRTModel_batch), 
                                                num_blocks = 4L)
                                                #num_blocks = ai(ceiling(N/mcmc_control$batch_size)))
        #Batch size controls computational cost and the variance of the log-likelihood estimate.
        #Block size (via num_blocks) controls the proposal variance for subsample updates (smaller blocks → gentler changes → higher acceptance rates).
      }
      if(mcmc_control$subsample_method == "full"){ 
        message("Enlisting NUTS kernels...")
        kernel <- lpme_env$numpyro$infer$NUTS(IRTModel_full, 
                                              max_tree_depth = ai(8), # 10 is default, slows performance down considerably
                                              target_accept_prob = 0.85 ) 
      }
      sampler <- lpme_env$numpyro$infer$MCMC(
        kernel,
        num_warmup = mcmc_control$n_samples_warmup,
        num_samples = mcmc_control$n_samples_mcmc,
        thinning = mcmc_control$n_thin_by, # Positive integer that controls the fraction of post-warmup samples that are retained. For example if thinning is 2 then every other sample is retained. Defaults to 1, i.e. no thinning.
        chain_method = mcmc_control$chain_method, # 'parallel' (default), 'sequential', 'vectorized'. 
        num_chains = mcmc_control$n_chains,
        jit_model_args = TRUE, 
        progress_bar = TRUE # set to TRUE for progress 
      )
      
      # run sampler with initialized abilities as COLMEANS of X (ASSUMPTION!)
      pdtype_ <- ddtype_ <- lpme_env$jnp$float32; lpme_env$jax$config$update("jax_enable_x64", FALSE) 
      #pdtype_ <- ddtype_ <- lpme_env$jnp$float64; lpme_env$jax$config$update("jax_enable_x64", TRUE) 
      #pdtype_ <- ddtype_ <- lpme_env$jnp$float16; lpme_env$jax$config$update("jax_enable_x64", FALSE) 
      
      # set initial parameters 
      UseEMInits <- FALSE
      ability_init <- lpme_env$jnp$broadcast_to(lpme_env$jnp$array(x_init<- (scale(rowMeans(observables_))*INIT_SCALER))$astype(pdtype_), list(mcmc_control$n_chains, N, 1L))
      if(UseEMInits){ ability_init <- lpme_env$jnp$broadcast_to(lpme_env$jnp$array( x_init<- scale(rowMeans(observables_))*INIT_SCALER )$astype(pdtype_), list(mcmc_control$n_chains, N, 1L)) } 
      
      difficulty_init <- lpme_env$jnp$array(matrix(rnorm(K*mcmc_control$n_chains,mean=0,sd=1/sqrt(K)),nrow=mcmc_control$n_chains) )$astype(pdtype_)
      if(UseEMInits){ difficulty_init <-  lpme_env$jnp$broadcast_to(lpme_env$jnp$array( out_emIRT$means$beta[,1]*INIT_SCALER )$astype(pdtype_), list(mcmc_control$n_chains, K)) }
      
      discrimination_init <- lpme_env$jnp$array(  matrix( (rnorm(K*mcmc_control$n_chains,mean=1,sd=1/sqrt(K))),nrow=mcmc_control$n_chains))$astype(pdtype_)
      if(UseEMInits){ discrimination_init <-  lpme_env$jnp$broadcast_to(lpme_env$jnp$array( out_emIRT$means$beta[,2] *INIT_SCALER )$astype(pdtype_),list(mcmc_control$n_chains, K)) }
      
      # run sampler 
      t0_ <- Sys.time()
      sampler$run(lpme_env$jax$random$PRNGKey( ai(runif(1,0,10000)) ), 
                  X = lpme_env$jnp$array(as.matrix(observables_))$astype( ddtype_ ),  # note: lpme_env$jnp$int16 here causes error (expects floats not ints)
                  Y = lpme_env$jnp$array(as.matrix(Y))$astype( ddtype_ ))
                  #init_params = list(
                  #  "ability" = ability_init,
                  #  "difficulty" = difficulty_init,
                  #  "discrimination" =  discrimination_init ) ) 
      PosteriorDraws <- sampler$get_samples(group_by_chain = TRUE)
      # PosteriorDraws$ability$shape; PosteriorDraws$ability[1,,1]
      ExtractAbil <- function(abil){ 
        abil <- do.call(cbind, sapply(1L:mcmc_control$n_chains, function(c_){
          abil_c <- as.matrix(lpme_env$np$array(abil)[c_,,,])
          return( list( t(abil_c) ) )
        }))
        theAnchor <- which.max(x_init)
        abil <- apply(abil,2,function(x){
                    if(x[theAnchor] < 0){ x <- -1*x } # anchor based on anchor 
                    return(x) })
        return( abil ) 
      }

      # summary( lm(Y~scale(AbilityMean) ) )
      # plot( as.matrix( lpme_env$np$array( PosteriorDraws$discrimination_oriented ) )[1,,1])
      # plot( as.matrix( lpme_env$np$array( PosteriorDraws$discrimination ) )[1,,1])
      # plot( as.matrix( lpme_env$np$array( PosteriorDraws$discrimination ) )[1,,2])
      
      # Calculate posterior means
      # ExtractAbil(PosteriorDraws$ability)
      # plot(ExtractAbil(PosteriorDraws$ability)[1,],main="i=1")
      # plot(ExtractAbil(PosteriorDraws$ability)[2,],main="i=1")
      AbilityMean <- rowMeans(  ExtractAbil(PosteriorDraws$ability) )
      DifficultyMean <- as.matrix(lpme_env$np$array(lpme_env$jnp$mean(PosteriorDraws$difficulty,0L:1L))) #  colMeans( as.matrix(lpme_env$np$array(PosteriorDraws$difficulty)) )
      
      # plot(scale(x.true[i_sampled]),scale(AbilityMean))
      # cor(x.true[i_sampled],AbilityMean)
      # plot(as.array(lpme_env$np$array(PosteriorDraws$ability))[1,721,])
      message(sprintf("\n MCMC Runtime: %.3f min",  tdiff_ <- as.numeric(difftime(Sys.time(),  t0_, units = "secs"))/60))
      message(sprintf("Mean(N-eff of nMCMC %%): %.2f%% \n", 100*mean(
              lpme_env$numpyro$diagnostics$effective_sample_size(# Computes effective sample size of input x, where the first dimension of x is chain dimension and the second dimension of x is draw dimension.
                lpme_env$jnp$reshape( PosteriorDraws$ability, list(mcmc_control$n_chains, ai(mcmc_control$n_samples_mcmc/mcmc_control$n_thin_by), N))
                ), na.rm=T)/(ai(mcmc_control$n_chains*mcmc_control$n_samples_mcmc/mcmc_control$n_thin_by) ) ))
      plot(lpme_env$np$array(PosteriorDraws$ability[1,,1,1]))
      if(estimation_method == "mcmc" & split_ == ""){ # method of compositions 
        # OuterNormed - this is what ideal does 
        RescaledAbilities_OuterNormed  <- (ExtractAbil(PosteriorDraws$ability)-mean(AbilityMean))/sd(AbilityMean)
        #Bayesian_OLSCoef_OuterNormed <- apply(RescaledAbilities, 2, function(x_){ coef(lm(Y~x_))[2]}) # simple MOC 
        Bayesian_OLSCoef_OuterNormed <- apply(RescaledAbilities_OuterNormed, 2, function(x_){ 
          VCovHat <- sandwich::vcovHC( myModel <- lm(Y~x_), type = "HC3" ) 
          coef_ <- mvtnorm::rmvnorm(n = 1, mean = coef(myModel), sigma = VCovHat)[1,2]
          } ) # complicated MOC 
        
        # summary( apply(RescaledAbilities_OuterNormed, 2, sd) ) 
        # sd(rowMeans(RescaledAbilities_OuterNormed));mean(rowMeans(RescaledAbilities_OuterNormed)) # confirm sanity value of 1, 0
        # hist( Bayesian_OLSCoef_OuterNormed ); summary( Bayesian_OLSCoef_OuterNormed)
        Bayesian_OLSSE_OuterNormed <- sd( Bayesian_OLSCoef_OuterNormed ) 
        Bayesian_OLSCoef_OuterNormed <- mean( Bayesian_OLSCoef_OuterNormed )
        
        # InnerNormed
        RescaledAbilities_InnerNormed  <- ( apply(ExtractAbil(PosteriorDraws$ability), 2, function(x_){scale(x_)})  ) 
        #Bayesian_OLSCoef_InnerNormed <- apply(RescaledAbilities_InnerNormed, 2, function(x_){ coef(lm(Y~x_))[2]}) # simple MOC 
        Bayesian_OLSCoef_InnerNormed <- apply(RescaledAbilities_InnerNormed, 2, function(x_){ 
          VCovHat <- sandwich::vcovHC( myModel <- lm(Y~x_), type = "HC3" ) 
          coef_ <- mvtnorm::rmvnorm(n = 1, mean = coef(myModel), sigma = VCovHat)[1,2]
        } ) # complicated MOC 
        # apply(RescaledAbilities_InnerNormed, 2, sd)
        # sd(rowMeans(RescaledAbilities_InnerNormed));mean(rowMeans(RescaledAbilities_InnerNormed)) # confirm sanity value of 1, 0
        # hist( Bayesian_OLSCoef_InnerNormed ); summary( Bayesian_OLSCoef_InnerNormed )
        Bayesian_OLSSE_InnerNormed <- sd( Bayesian_OLSCoef_InnerNormed ) 
        Bayesian_OLSCoef_InnerNormed <- mean( Bayesian_OLSCoef_InnerNormed )
      }
      if(estimation_method == "mcmc_joint" & split_ == ""){ # full bayesian model 
        #RescaledAbilities <- (ExtractAbil(PosteriorDraws$ability)-mean(AbilityMean))/sd(AbilityMean)

        # note: multiplication of coefficient by sd generates interpretation of coeff 
        # as representing change in outcome associated with 1 unit deviation
        Bayesian_OLSCoef_OuterNormed <- c(as.matrix(lpme_env$np$array(PosteriorDraws$YModel_slope))) * sd(AbilityMean)
        # sd(rowMeans(RescaledAbilities));mean(rowMeans(RescaledAbilities)) # confirm sanity values of 1,0 
        # hist(Bayesian_OLSCoef_OuterNormed);abline(v=0.4); summary( Bayesian_OLSCoef_OuterNormed )
        #Bayesian_OLSCoef_OuterNormed[Bayesian_OLSCoef_OuterNormed<0] <- -1*Bayesian_OLSCoef_OuterNormed[Bayesian_OLSCoef_OuterNormed<0]
        Bayesian_OLSSE_OuterNormed <- sd( Bayesian_OLSCoef_OuterNormed ) 
        Bayesian_OLSCoef_OuterNormed <- mean( Bayesian_OLSCoef_OuterNormed )
        
        InnerSDs <- apply(ExtractAbil(PosteriorDraws$ability), 2, sd)
        Bayesian_OLSCoef_InnerNormed <- c(as.matrix(lpme_env$np$array(PosteriorDraws$YModel_slope))) * InnerSDs
        # sd(rowMeans(RescaledAbilities));mean(rowMeans(RescaledAbilities)) # confirm sanity values of 1,0 
        # hist(Bayesian_OLSCoef_InnerNormed); abline(v=0.4,lwd=2); summary( Bayesian_OLSCoef_InnerNormed )
        Bayesian_OLSSE_InnerNormed <- sd( Bayesian_OLSCoef_InnerNormed ) 
        Bayesian_OLSCoef_InnerNormed <- mean( Bayesian_OLSCoef_InnerNormed )
      }
      
      # rescale 
      x.est_MCMC <- x.est_ <- as.matrix(scale(AbilityMean)); s_past <- 1
      if(estimation_method == "mcmc_overimputation" & split_ == ""){ 
        for(outType_ in c("Outer","Inner")){ 
          # file:///Users/cjerzak/Dropbox/LatentMeasures/literature/CAUGHEY-ps8-solution.html
          if(outType_ == "Outer"){
            RescaledAbilities  <- (ExtractAbil(PosteriorDraws$ability)-mean(AbilityMean))/sd(AbilityMean)
          }
          if(outType_ == "Inner"){
            RescaledAbilities  <- ( apply(ExtractAbil(PosteriorDraws$ability), 2, function(x_){scale(x_)})  ) 
          }
          Xobs_mean <- apply(RescaledAbilities, 1, function(x_){ mean(x_) } ) 
          Xobs_SE <- apply(RescaledAbilities, 1, function(x_){ sd(x_) } ) 
          
          dat_ <- cbind(Y, x.est_MCMC)
          
          # Specify (over-)imputation model
          # priors: #a numeric matrix with four columns 
          # (row, column, mean, standard deviation) 
          # indicating noisy estimates of the values to be imputed.
          outcome_priors <- cbind(
            1:nrow(dat_), 1,              
            Y, # mean 
            1 # sd 
          )
          policy_priors <- cbind( 
            1:nrow(dat_), 2,              
            Xobs_mean,
            Xobs_SE
          )
          #priors <- rbind(policy_priors, outcome_priors)
          priors <- policy_priors
          
          #a numeric matrix where each row indicates a row and column of x to be overimputed.
          # overimp <- rbind(cbind(1:nrow(dat_), 1), cbind(1:nrow(dat_), 2))
          overimp <- cbind(1:nrow(dat_), 2)
          
          # perform overimputation 
          overimputed_data <- Amelia::amelia( 
            x = dat_, m = (nOverImpute <- 5L),
            p2s = 0,
            priors = priors, 
            overimp = overimp 
          )
          
          # Perform multiple overimputation
          overimputed_Y <- do.call(cbind,lapply(overimputed_data$imputations,function(l_){l_[,1]}))
          overimputed_x.est_MCMC <- do.call(cbind,lapply(overimputed_data$imputations,function(l_){l_[,2]}))
          if(outType_ == "Outer"){
            overimputed_x.est_MCMC  <- (overimputed_x.est_MCMC-mean(rowMeans(overimputed_x.est_MCMC)))/sd(rowMeans(overimputed_x.est_MCMC))
            #sd(rowMeans(overimputed_x.est_MCMC))
          }
          if(outType_ == "Inner"){
            overimputed_x.est_MCMC  <- ( apply(overimputed_x.est_MCMC, 2, function(x_){scale(x_)})  ) 
            #apply(overimputed_x.est_MCMC,2,sd)
          }
          
          # cor(cbind(x.est_MCMC, overimputed_x.est_MCMC))
          # plot(x.est_MCMC,overimputed_x.est_MCMC[,1])
          
          # Analyze overimputed datasets
          overimputed_coefs <- unlist(unlist( sapply(1:nOverImpute, function(s_){
            coef(lm(overimputed_Y[,s_] ~ overimputed_x.est_MCMC[,s_]))[2]
          })))
          if(outType_ == "Outer"){ 
            Bayesian_OLSSE_OuterNormed <- sd( overimputed_coefs )
            Bayesian_OLSCoef_OuterNormed <- mean( overimputed_coefs )
          }
          if(outType_ == "Inner"){ 
            Bayesian_OLSSE_InnerNormed <- sd( overimputed_coefs )
            Bayesian_OLSCoef_InnerNormed <- mean( overimputed_coefs )
          }
        }
      }
    }
    
    if(split_ %in% c("1","2")){
      if(cor(x.est_, x_est_past) < 0){
        x.est_  <- -1*x.est_
      }
    }
    x_est_past <- x.est_
    eval(parse(text = sprintf("x.est%s <- x.est_", split_)))
  }
  # c(mean(x.est),sd(x.est))
  
  # simple linear reg 
  theOLS <- lm(Y ~ x.est)
  
  # stage 1 results
  IVStage1 <- lm(x.est2 ~ x.est1)
  
  # corrected IV
  # cor(x.est1,x.est2)
  IVStage2_a <- AER::ivreg(Y ~ x.est2 | x.est1)
  IVStage2_b <- AER::ivreg(Y ~ x.est1 | x.est2)
  IVRegCoef_a <- coef(IVStage2_a)[2]
  IVRegCoef_b <- coef(IVStage2_b)[2]
  IVRegCoef <- (IVRegCoef_a + IVRegCoef_b) / 2
  Corrected_IVRegCoef_a <- (coef(IVStage2_a)[2] * sqrt( max(c((ep_ <- 0.01), 
                                                              cor(x.est1, x.est2) ))) )
  Corrected_IVRegCoef_b <- (coef(IVStage2_b)[2] * sqrt( max(c(ep_,
                                                              cor(x.est1, x.est2) ))) )
  Corrected_IVRegCoef <- ( Corrected_IVRegCoef_a + Corrected_IVRegCoef_b )/2
  
  # corrected OLS 
  Corrected_OLSCoef_a <- coef(lm(Y ~ x.est1))[2] * (CorrectionFactor <- 1/sqrt( max(c(ep_, cor(x.est1, x.est2) ))))
  Corrected_OLSCoef_b <- coef(lm(Y ~ x.est2))[2] * CorrectionFactor
  Corrected_OLSCoef <- (Corrected_OLSCoef_a + Corrected_OLSCoef_b)/2
  t_OneRun <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  # ERV analysis 
  mstage1 <- lm(x.est2 ~ x.est1)
  mreduced <- lm(Y ~ x.est1)
  mstage1ERV <- sensemakr::extreme_robustness_value( mstage1 )[[2]]
  mreducedERV <- sensemakr::extreme_robustness_value( mreduced )[[2]]
  
  # save results 
  results <- list(
    "ols_coef" = coef(summary(theOLS))[2, 1],
    "ols_se" = coef(summary(theOLS))[2, 2],
    "ols_tstat" = coef(summary(theOLS))[2, 3],
    
    "iv_coef_a" = IVRegCoef_a, 
    "iv_coef_b" = IVRegCoef_b, 
    "iv_coef" = IVRegCoef, 
    "iv_se" = coef(summary(IVStage2_a))[2, 2],
    "iv_tstat" = coef(summary(IVStage2_a))[2, 3],
    
    "corrected_iv_coef_a" = Corrected_IVRegCoef_a,
    "corrected_iv_coef_b" = Corrected_IVRegCoef_b,
    "corrected_iv_coef" = Corrected_IVRegCoef,
    "corrected_iv_se" = NA,
    "corrected_iv_tstat" = Corrected_IVRegCoef / coef(summary(IVStage2_a))[2, 2],
    "var_est_split" = var(x.est1 - x.est2) / 2,  # var_est_split
    
    "corrected_ols_coef_a" = Corrected_OLSCoef_a, 
    "corrected_ols_coef_b" = Corrected_OLSCoef_b, 
    "corrected_ols_coef" = Corrected_OLSCoef, 
    "corrected_ols_se" = NA,
    "corrected_ols_tstat" = NA,
    
    "corrected_ols_coef_alt" = NA, 
    "corrected_ols_se_alt" = NA,
    "corrected_ols_tstat_alt" = NA,
    
    "bayesian_ols_coef_outer_normed" = Bayesian_OLSCoef_OuterNormed, 
    "bayesian_ols_se_outer_normed" = Bayesian_OLSSE_OuterNormed,
    "bayesian_ols_tstat_outer_normed" = Bayesian_OLSCoef_OuterNormed / Bayesian_OLSSE_OuterNormed,
    
    "bayesian_ols_coef_inner_normed" = Bayesian_OLSCoef_InnerNormed, 
    "bayesian_ols_se_inner_normed" = Bayesian_OLSSE_InnerNormed,
    "bayesian_ols_tstat_inner_normed" = Bayesian_OLSCoef_InnerNormed / Bayesian_OLSSE_InnerNormed,
    
    "m_stage_1_erv" = mstage1ERV, 
    "m_reduced_erv" = mreducedERV, 
    
    "x_est1" = x.est1,
    "x_est2" = x.est2
  )
    class(results) <- "lpme_onerun"
    return( results ) 
}
