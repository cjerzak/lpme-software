#!/usr/bin/env Rscript
#' lpme_OneRun
#'
#' Implements analysis for latent variable models with measurement error correction
#'
#' @param Yobs A vector of observed outcome variables
#' @param ObservablesMat A matrix of observable indicators used to estimate the latent variable
#' @param ObservablesGroupings A vector specifying groupings for the observable indicators. Default is column names of ObservablesMat.
#' @param MakeObservablesGroupings Logical. If TRUE, creates dummy variables for each level of the observable indicators. Default is FALSE.
#' @param seed Random seed for reproducibility. Default is a random integer between 1 and 10000.
#'
#' @return A list containing various estimates and statistics:
#' \itemize{
#'   \item OLSCoef: Coefficient from naive OLS regression
#'   \item OLSSE: Standard error of naive OLS coefficient
#'   \item OLSTstat: T-statistic of naive OLS coefficient
#'   \item Corrected_OLSCoef: OLS coefficient corrected for measurement error
#'   \item Corrected_OLSSE: Standard error of corrected OLS coefficient (currently NA)
#'   \item Corrected_OLSTstat: T-statistic of corrected OLS coefficient (currently NA)
#'   \item Corrected_OLSCoef_alt: Alternative corrected OLS coefficient
#'   \item IVRegCoef: Coefficient from instrumental variable regression
#'   \item IVRegSE: Standard error of IV regression coefficient
#'   \item IVRegTstat: T-statistic of IV regression coefficient
#'   \item x.est1: First set of latent variable estimates
#'   \item x.est2: Second set of latent variable estimates
#'   \item Corrected_IVRegCoef: IV regression coefficient corrected for measurement error
#'   \item Corrected_IVRegSE: Standard error of corrected IV coefficient (currently NA)
#'   \item Corrected_IVRegTstat: T-statistic of corrected IV coefficient
#'   \item VarEst_split: Estimated variance of the measurement error
#' }
#'
#' @details 
#' This function implements a latent variable analysis with measurement error correction. 
#' It splits the observable indicators into two sets, estimates latent variables using each set, 
#' and then applies various correction methods including OLS correction and instrumental variable approaches.
#'
#' @examples
#' # Generate some example data
#' set.seed(123)
#' Yobs <- rnorm(100)
#' ObservablesMat <- as.data.frame( matrix(sample(c(0,1), 1000*10, replace = T), ncol = 10) )
#' 
#' # Run the analysis
#' results <- lpme_OneRun(Yobs, ObservablesMat)
#' 
#' # View the corrected OLS coefficient
#' print(results$Corrected_OLSCoef)
#'
#' @export
#' @importFrom stats lm cor var rnorm
#' @importFrom AER ivreg
#' @importFrom MCMCpack rollcall

lpme_OneRun <- function(Yobs,
                         ObservablesMat, 
                         ObservablesGroupings = colnames(ObservablesMat),
                         MakeObservablesGroupings = F, 
                         EstimationMethod = "emIRT", 
                         conda_env = NULL, 
                         Sys.setenv_text = NULL,
                         seed = NULL){
  library(emIRT);
  starting_seed <- sample(runif(1,1,10000))
  FullBayesiasn_OLSCoef_se <- FullBayesiasn_OLSCoef <- NA; 
  if(!is.null(seed)){ set.seed(seed) } 
  items.split1_names <- sample(unique(ObservablesGroupings), 
                               size = floor(length(unique(ObservablesGroupings))/2), replace=F)
  items.split2_names <- unique(ObservablesGroupings)[! (ObservablesGroupings %in% items.split1_names)]
  for(split_ in c("", "1", "2")){
    if(split_ == ""){ items.split_ <- 1:length(ObservablesGroupings) }
    if(split_ == "1"){ items.split_ <- (1:length(ObservablesGroupings))[ObservablesGroupings %in% items.split1_names] }
    if(split_ == "2"){ items.split_ <- (1:length(ObservablesGroupings))[ObservablesGroupings %in% items.split2_names] }
    
    # estimating ideal points
    if(MakeObservablesGroupings == F){
      ObservablesMat_ <- ObservablesMat[,items.split_]
    }
    if(MakeObservablesGroupings == T){
      ObservablesMat_ <- do.call(cbind,unlist(apply(ObservablesMat[,items.split_],2,function(zer){
        list( model.matrix(~0+as.factor(zer))[,-1] )}),recursive = F))
    }
    
    if(grepl(EstimationMethod,pattern = "MCMC")){
      if(split_ == ""){
        # Load required libraries
        library(reticulate)
        reticulate::use_condaenv(conda_env)
        
        # set environmental variables 
        if( !is.null(Sys.setenv_text) ){ 
          eval(parse(text = Sys.setenv_text), envir = .GlobalEnv)
        }
  
        # Import necessary Python modules
        np <- import("numpy", convert = F)
        jax <- import("jax")
        jnp <- import("jax.numpy")
        random <- import("jax.random")
        numpyro <- import("numpyro")
        dist <- import("numpyro.distributions")
        f2i <- function(f_){jnp$array(f_,jnp$int32)}
        f2a <- function(x){jnp$array(x,jnp$float32)}
        ai <- as.integer
    
    }  
        # Construct for annotating conditionally independent variables.
        # Within a plate context manager, sample sites will be automatically broadcasted
        # to the size of the plate. 
        # Additionally, a scale factor might be applied by certain inference algorithms
        # if  subsample_size is specified.
        
        # Set up MCMC
        nChains <- 2L
        numpyro$set_host_device_count(nDevices <- 1L)
        ChainMethod <- "sequential"
        nSamplesWarmup <- (1000L)
        nSamplesMCMC <- (500L)
        nThinBy <- 2L
        N <- ai(nrow(ObservablesMat_))
        K <- ai(ncol(ObservablesMat_))
        
        # EstimationMethod <- "MCMCFull"
        # EstimationMethod <- "MCMC"
        
        # Define the two-parameter IRT model
        irt_model <- function(X, # binary indicators 
                              Y, # outcome (used if EstimationMethod <- "MCMCFull")
                              N, # number of observations  
                              K # number of items 
                              ) {
          # Priors
          with(numpyro$plate("rows", N), {
            ability <- numpyro$sample("ability", dist$Normal(0, 1))
          })
          with(numpyro$plate("columns", K), {
            difficulty <- numpyro$sample("difficulty", dist$Normal(0, 1))
            discrimination <- numpyro$sample("discrimination", dist$LogNormal(0, 1)) # mass on positive values 
          })
  
          # likelihood
          logits <- jnp$outer(ability, discrimination) - difficulty
          
          # sanity check prints 
          if(T == F){ 
            print("ability shape:"); print(ability$shape)
            print("discrimination shape:"); print(discrimination$shape)
            print("difficulty shape:");print(difficulty$shape)
            print("X shape:");print(X$shape)
            print("logits shape:");print(logits$shape)
          }
          
          numpyro$sample("Xlik", dist$Bernoulli(logits = logits), obs = X)
          
          if(EstimationMethod == "MCMCFull"){
            # Outcome priors 
            Y_intercept <- numpyro$sample("YModel_intercept", dist$Normal(0, 1))
            #Y_slope <- numpyro$sample("YModel_slope", dist$Normal(0, 2))
            Y_slope <- numpyro$sample("YModel_slope", dist$LogNormal(0, 1)) # constrain slope to be positive 
            #Y_sigma <- numpyro$sample("YModel_sigma", dist$LogNormal(0, 1))
            Y_sigma <- numpyro$sample("YModel_sigma", dist$HalfNormal(1))
            
            # Outcome model likelihood
            ability <- jnp$expand_dims(ability,1L)
            Y_mu <- Y_intercept + jnp$multiply(Y_slope, jax$nn$standardize(ability,0L) )
            
            if(T == F){ # sanity check prints 
              print("ability");print(ability$shape)
              print("Y_mu"); print(Y_mu$shape)
              print("Y"); print(Y$shape)
              print("Y_slope"); print(Y_slope$shape)
              print("Y_sigma"); print(Y_sigma$shape)
            }
            numpyro$sample("Ylik", dist$Normal(Y_mu, Y_sigma), obs=Y) # Y_mu has to have same shape as Y
          }
        }
        
      # setup & run MCMC 
      sampler <- numpyro$infer$MCMC(
        numpyro$infer$NUTS(irt_model),
        num_warmup = nSamplesWarmup,
        num_samples = nSamplesMCMC,
        thinning = nThinBy, # Positive integer that controls the fraction of post-warmup samples that are retained. For example if thinning is 2 then every other sample is retained. Defaults to 1, i.e. no thinning.
        chain_method = ChainMethod, # ‘parallel’ (default), ‘sequential’, ‘vectorized’. 
        num_chains = nChains
      )
      sampler$run(jax$random$PRNGKey( ai(runif(1,0,10000)) ), 
                  X = jnp$array(as.matrix(ObservablesMat_))$astype(jnp$int16), 
                  Y = jnp$array(as.matrix(Yobs))$astype(jnp$float32), 
                  N = N, K = K) # don't use numpy array for N and K inputs here! 
      PosteriorDraws <- sampler$get_samples(group_by_chain = T)
      if(EstimationMethod == "MCMCFull" & split_ == ""){
        FullBayesiasn_OLSCoef <- c(as.matrix(np$array(PosteriorDraws$YModel_slope)))
        # hist(FullBayesiasn_OLSCoef)
        FullBayesiasn_OLSSE <- sd( FullBayesiasn_OLSCoef )
        FullBayesiasn_OLSCoef <- mean( FullBayesiasn_OLSCoef )
      }
      
      ExtractAbil <- function(abil){ # note: deals with identifiability of scale 
        abil <- do.call(cbind, sapply(0L:(nChains-1L), function(c_){
          abil_c <- as.matrix(np$array(abil)[c_,,])
          abil_c <- t(abil_c)
          if(cor(rowMeans(abil_c),Yobs)<0){ abil_c <- -1*abil_c }
          return( list(abil_c) )
        }))
        return( abil ) 
      }
      
      # Calculate posterior means
      AbilityMean <- rowMeans(  ExtractAbil(PosteriorDraws$ability) )
      DifficultyMean <- as.matrix(np$array(jnp$mean(PosteriorDraws$difficulty,0L:1L))) #  colMeans( as.matrix(np$array(PosteriorDraws$difficulty)) )
      DiscriminationMean <- as.matrix(np$array(jnp$mean(PosteriorDraws$discrimination,0L:1L))) # colMeans( as.matrix(np$array(PosteriorDraws$discrimination)) )
      
      # resacle 
      x.est_MCMC <- x.est_ <- as.matrix(scale(AbilityMean)); s_past <- 1
      # as.matrix(np$array(PosteriorDraws$ability))[,1:3]
      # plot(DifficultyMean,DiscriminationMean)
      # plot(x.est_MCMC, x.est_EM); abline(a=0,b=1); cor(x.est_MCMC, x.est_EM)
    }
      
    if(EstimationMethod == "emIRT"){ 
      x_init <- apply( ObservablesMat_, 1, function(x){ mean(f2n(x), na.rm=T)})
      rc_ <- convertRC( rollcall(ObservablesMat_) )
      #s_ <- list("alpha" = matrix(rnorm(ncol(ObservablesMat_), sd = 1)),"beta" = matrix(rnorm(ncol(ObservablesMat_), sd = 1)), "x" = matrix(x_init))
      s_ <- getStarts(.N= rc_$n, .J = rc_$m, .D = 1)
      
      # fixing in case directions of s1 and s2 x starts are flipped
      if(split_ %in% c("1","2")){ if(cor(s_$x, s_past$x) < 0){ s_$x <- -s_$x; s_$beta <- -s_$beta } }
      lout.sim_ <- binIRT(.rc = rc_, 
                          .starts = s_, 
                          .priors = makePriors(.N= rc_$n, .J = rc_$m, .D = 1), 
                          .control= list(threads=1, verbose=FALSE, thresh=1e-6,verbose=F),
                          .anchor_subject = which.max(x_init)) # set direction
      x.est_EM <- x.est_ <- scale(lout.sim_$means$x); s_past <- s_
    }
    if(cor(x.est_, Yobs, use="p") < 0){ x.est_ <- -x.est_ }
    eval(parse(text = sprintf("x.est%s <- x.est_", split_)))
  }
  if(!is.null(seed)){ set.seed(starting_seed) } 

  # simple linear reg 
  theOLS <- lm(Yobs ~ x.est)
  
  # stage 1 results
  IVStage1 <- lm(x.est2 ~ x.est1)
  
  # cor epsilon
  ep_ <- 0.01
  
  # corrected IV
  IVStage2_a <- AER::ivreg(Yobs ~ x.est2 | x.est1)
  IVStage2_b <- AER::ivreg(Yobs ~ x.est1 | x.est2)
  Corrected_IVRegCoef_a <- (coef(IVStage2_a)[2] * sqrt( max(c(ep_, cor(x.est1, x.est2) ))))
  Corrected_IVRegCoef_b <- (coef(IVStage2_b)[2] * sqrt( max(c(ep_, cor(x.est1, x.est2) ))))
  Corrected_IVRegCoef <- ( Corrected_IVRegCoef_a + Corrected_IVRegCoef_b )/2
  
  # corrected OLS 
  Corrected_OLSCoef1 <- coef(lm(Yobs ~ x.est1))[2] * (CorrectionFactor <- 1/sqrt( max(c(ep_, cor(x.est1, x.est2) ))))
  Corrected_OLSCoef2 <- coef(lm(Yobs ~ x.est2))[2] * CorrectionFactor
  Corrected_OLSCoef <- (Corrected_OLSCoef1 + Corrected_OLSCoef2)/2
  # 0.3923193 
  
  # ERV analysis 
  mstage1 <- lm(x.est2 ~ x.est1)
  mreduced <- lm(Yobs ~ x.est1)
  mstage1ERV <- sensemakr::extreme_robustness_value( mstage1 )[[2]]
  mreducedERV <- sensemakr::extreme_robustness_value( mreduced )[[2]]
  
  # save results 
  ret_ <- list("OLSCoef" = coef(summary(theOLS))[2,1],
       "OLSSE" = coef(summary(theOLS))[2,2],
       "OLSTstat" = coef(summary(theOLS))[2,3],
       
       "IVRegCoef" = coef(summary(IVStage2_a))[2,1], 
       "IVRegSE" = coef(summary(IVStage2_a))[2,2],
       "IVRegTstat" = coef(summary(IVStage2_a))[2,3],
       
       "Corrected_IVRegCoef" = Corrected_IVRegCoef,
       "Corrected_IVRegSE" = NA,
       "Corrected_IVRegTstat" =  Corrected_IVRegCoef / coef(summary(IVStage2_a))[2,2],
       "VarEst_split" = var(x.est1 - x.est2) / 2 ,
  
        "Corrected_OLSCoef" = Corrected_OLSCoef, 
        "Corrected_OLSSE" = NA,
        "Corrected_OLSTstat" = NA,
        
        "Corrected_OLSCoef_alt" = NA, 
        "Corrected_OLSSE_alt" = NA,
        "Corrected_OLSTstat_alt" = NA,
       
       "FullBayesiasn_OLSCoef" = FullBayesiasn_OLSCoef, 
       "FullBayesiasn_OLSSE" = FullBayesiasn_OLSSE,
       "FullBayesiasn_OLSTstat" = FullBayesiasn_OLSCoef/FullBayesiasn_OLSSE,
    
        "mstage1ERV" = mstage1ERV, 
        "mreducedERV" = mreducedERV, 
  
        "x.est1" = x.est1,
        "x.est2" = x.est2)

    return( ret_ ) 
}
