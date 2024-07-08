#!/usr/bin/env Rscript
#' lpme
#'
#' Implements bootstrapped analysis for latent variable models with measurement error correction
#'
#' @param Yobs A vector of observed outcome variables
#' @param ObservablesMat A matrix of observable indicators used to estimate the latent variable
#' @param ObservablesGroupings A vector specifying groupings for the observable indicators. Default is column names of ObservablesMat.
#' @param MakeObservablesGroupings Logical. If TRUE, creates dummy variables for each level of the observable indicators. Default is FALSE.
#' @param nBoot Integer. Number of bootstrap iterations. Default is 32.
#' @param nPartition Integer. Number of partitions for each bootstrap iteration. Default is 10.
#' @param bootBasis Vector of indices or grouping variable for stratified bootstrap. Default is 1:length(Yobs).
#' @param ReturnIntermediaries Logical. If TRUE, returns intermediate results. Default is TRUE.
#' @param seed Random seed for reproducibility. Default is a random integer between 1 and 10000.
#'
#' @return A list containing various estimates and statistics:
#' \itemize{
#'   \item OLSCoef: Coefficient from naive OLS regression
#'   \item OLSSE: Bootstrap standard error of naive OLS coefficient
#'   \item OLSTstat: T-statistic of naive OLS coefficient
#'   \item Corrected_OLSCoef: OLS coefficient corrected for measurement error
#'   \item Corrected_OLSSE: Bootstrap standard error of corrected OLS coefficient
#'   \item Corrected_OLSTstat: T-statistic of corrected OLS coefficient
#'   \item Corrected_OLSCoef_alt: Alternative corrected OLS coefficient
#'   \item Corrected_OLSSE_alt: Bootstrap standard error of alternative corrected OLS coefficient
#'   \item Corrected_OLSTstat_alt: T-statistic of alternative corrected OLS coefficient
#'   \item IVRegCoef: Coefficient from instrumental variable regression
#'   \item IVRegSE: Bootstrap standard error of IV regression coefficient
#'   \item IVRegTstat: T-statistic of IV regression coefficient
#'   \item Corrected_IVRegCoef: IV regression coefficient corrected for measurement error
#'   \item Corrected_IVRegSE: Bootstrap standard error of corrected IV coefficient
#'   \item Corrected_IVRegTstat: T-statistic of corrected IV coefficient
#'   \item x.est1: First set of latent variable estimates
#'   \item x.est2: Second set of latent variable estimates
#'   \item VarEst_split: Estimated variance of the measurement error
#'   \item VarEst_split_se: Bootstrap standard error of the estimated variance of the measurement error
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
#' Yobs <- rnorm(1000)
#' ObservablesMat <- as.data.frame( matrix(sample(c(0,1), 1000*10, replace = T), ncol = 10) )
#' 
#' # Run the bootstrapped analysis
#' results <- lpme(Yobs, ObservablesMat, nBoot = 100, nPartition = 5)
#' 
#' # View the corrected IV coefficient and its standard error
#' print(c(results$Corrected_IVRegCoef, results$Corrected_IVRegSE))
#'
#' @export
#' @importFrom stats sd median tapply
#' @importFrom lpme lpme_OneRun

lpme <- function(Yobs,
                      ObservablesMat, 
                      ObservablesGroupings = colnames(ObservablesMat),
                      MakeObservablesGroupings = F,
                      nBoot = 32L, nPartition = 10L, 
                      bootBasis = 1:length(Yobs),
                      ReturnIntermediaries = T, 
                      seed = runif(1, 1, 10000)){ 
  set.seed(seed)
  for(booti_ in 1L:(nBoot+1L)){
    boot_indices <- 1:length(Yobs); if(booti_ > 1L){
        boot_indices <- sample(unique(as.character(bootBasis)),length(unique(bootBasis)), replace = T)
        boot_indices <- unlist(tapply(1:length(bootBasis),as.character(bootBasis),c)[boot_indices])
    }
    for(parti_ in 1L:nPartition){
      LatentRunResults_ <- lpme::lpme_OneRun(
                   Yobs[boot_indices],
                   ObservablesMat[boot_indices,], 
                   ObservablesGroupings = colnames(ObservablesMat),
                   MakeObservablesGroupings = F , 
                   seed = seed+parti_)
      
      # save indices indices 
      LatentRunResults_$PartitionIndex <- parti_; LatentRunResults_$BootIndex <- booti_
      
      # store results 
      if(booti_ == 1 & parti_ == 1){  LatentRunResults <- LatentRunResults_  }
      if(!(booti_ == 1 & parti_ == 1)){ for(name_ in names(LatentRunResults_)){ 
        eval(parse(text = sprintf("LatentRunResults$%s <- cbind(LatentRunResults$%s, LatentRunResults_$%s)",name_,name_,name_ )))
      } }
    } }
  theSumFxn <- median 
  #theSumFxn <- mean
  names( LatentRunResults ) <- paste0("Intermediary_",names(LatentRunResults))
  VarEst_split <- theSumFxn(apply(LatentRunResults$Intermediary_x.est1[,which(LatentRunResults$Intermediary_BootIndex==1)] - 
                             LatentRunResults$Intermediary_x.est2[,which(LatentRunResults$Intermediary_BootIndex==1)], 1, sd))
  VarEst_split_se <- sd( sapply(2:(nBoot+1),function(boot_){
                            theSumFxn(apply(LatentRunResults$Intermediary_x.est1[,which(LatentRunResults$Intermediary_BootIndex==boot_)] - 
                               LatentRunResults$Intermediary_x.est2[,which(LatentRunResults$Intermediary_BootIndex==boot_)], 1, sd)) }))
    
  return( 
    list(
       "OLSCoef" = (m1_ <- tapply(LatentRunResults$Intermediary_OLSCoef,LatentRunResults$Intermediary_BootIndex,theSumFxn)[1]),
       "OLSSE" = (se1_ <- sd( tapply(LatentRunResults$Intermediary_OLSCoef,LatentRunResults$Intermediary_BootIndex,theSumFxn)[-1] )),
       "OLSTstat" = (m1_/se1_),
       
       "Corrected_OLSCoef" = (m1b_ <- tapply(LatentRunResults$Intermediary_Corrected_OLSCoef,LatentRunResults$Intermediary_BootIndex,theSumFxn)[1]),
       "Corrected_OLSSE" = (se1b_ <- sd( tapply(LatentRunResults$Intermediary_Corrected_OLSCoef,LatentRunResults$Intermediary_BootIndex,theSumFxn)[-1] )),
       "Corrected_OLSTstat" = (m1b_/se1b_),
       
       "Corrected_OLSCoef_alt" = (m1b_ <- tapply(LatentRunResults$Intermediary_Corrected_OLSCoef_alt,LatentRunResults$Intermediary_BootIndex,theSumFxn)[1]),
       "Corrected_OLSSE_alt" = (se1b_ <- sd( tapply(LatentRunResults$Intermediary_Corrected_OLSCoef_alt,LatentRunResults$Intermediary_BootIndex,theSumFxn)[-1] )),
       "Corrected_OLSTstat_alt" = (m1b_/se1b_),
       
       "IVRegCoef" = (m2_ <- tapply(LatentRunResults$Intermediary_IVRegCoef,LatentRunResults$Intermediary_BootIndex,theSumFxn)[1]),
       "IVRegSE" = (se2_ <- sd(tapply(LatentRunResults$Intermediary_IVRegCoef,LatentRunResults$Intermediary_BootIndex,theSumFxn)[-1] )),
       "IVRegTstat" = (m2_/se2_),
       
       "Corrected_IVRegCoef" = (m4_ <- tapply(LatentRunResults$Intermediary_Corrected_IVRegCoef,LatentRunResults$Intermediary_BootIndex,theSumFxn)[1]),
       "Corrected_IVRegSE" = (se4_ <- sd(tapply(LatentRunResults$Intermediary_Corrected_IVRegCoef,LatentRunResults$Intermediary_BootIndex,theSumFxn)[-1] )),
       "Corrected_IVRegTstat"  = (m4_/se4_),
       
       "x.est1" = LatentRunResults$Intermediary_x.est1[,1],
       "x.est2" = LatentRunResults$Intermediary_x.est2[,1],
       
       "VarEst_split" = VarEst_split,
       "VarEst_split_se" = VarEst_split_se)
  ) 
}
