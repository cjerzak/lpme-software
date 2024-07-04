#!/usr/bin/env Rscript
#' LatentErrorBoot
#'
#' Implements fancy analysis
#'
#' @param x The input
#' @return `x` Returns this. 
#' @export
#'
#' @details Details. 
#'
#' @examples
#'
#' # Comment here
#' #x <- LatentErrorBoot(x)
#' @export
#' @md

LatentRun <- function(Yobs,
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
      LatentRunResults_ <- latenterror::LatentOneRun(
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

# core estimation function base model
# estimation function
# -- how many approaches to fix 
# -- estimate degree of measurement error 
# -- how much does split half affect us? 
# -- results about how things differ in latent case 
# correction functions 
# -- 
# diagnostic function

