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

#' @import Rfast
#' @md

LatentRun <- function(Yobs,
                      ObservablesMat, 
                      ObservablesGroupings = colnames(ObservablesMat),
                      MakeObservablesGroupings = F,
                      nBoot = 32L, nPartition = 10L, 
                      ReturnIntermediaries = T){ 
  for(booti_ in 1L:(nBoot+1L)){
    for(parti_ in 1L:nPartition){
      boot_indices <- 1:length(Yobs); if(booti_ > 1L){
        boot_indices <- sample(1:length(Yobs),length(Yobs),replace=T)
      }
      LatentRunResults_ <- latenterror::LatentOneRun(Yobs[boot_indices],
                   ObservablesMat[boot_indices,], 
                   ObservablesGroupings = colnames(ObservablesMat),
                   MakeObservablesGroupings = F)
      
      # save indices indices 
      LatentRunResults_$PartitionIndex <- parti_; LatentRunResults_$BootIndex <- booti_
      
      # store results 
      if(booti_ == 1 & parti_ == 1){  LatentRunResults <- LatentRunResults_  }
      if(!(booti_ == 1 & parti_ == 1)){ for(name_ in names(LatentRunResults_)){ 
        eval(parse(text = sprintf("LatentRunResults$%s <- cbind(LatentRunResults$%s,
                                                                LatentRunResults_$%s)",name_,name_,name_ )))
      } }
    } }
  names( LatentRunResults ) <- paste0("Intermediary_",names(LatentRunResults))
  VarEst_split <- mean(apply(LatentRunResults$Intermediary_x.est1[,which(LatentRunResults$Intermediary_BootIndex==1)] - 
                             LatentRunResults$Intermediary_x.est2[,which(LatentRunResults$Intermediary_BootIndex==1)], 1, sd))
  VarEst_split_se <- sd( sapply(2:(nBoot+1),function(boot_){
                            mean(apply(LatentRunResults$Intermediary_x.est1[,which(LatentRunResults$Intermediary_BootIndex==boot_)] - 
                               LatentRunResults$Intermediary_x.est2[,which(LatentRunResults$Intermediary_BootIndex==boot_)], 1, sd)) }))
  Corrected_IVRegCoef_vec <- tapply( (LatentRunResults$Intermediary_IVRegCoef / sqrt( 1 + var( 
                  apply(LatentRunResults$Intermediary_x.est1 - LatentRunResults$Intermediary_x.est2,1,sd) )/2 )),  
                  LatentRunResults$Intermediary_BootIndex, mean)
    
  LatentRunResults <- list("OLSCoef" = (m1_ <- tapply(LatentRunResults$Intermediary_OLSCoef,LatentRunResults$Intermediary_BootIndex,mean)[1]),
       "OLSSE" = (se1_ <- sd( tapply(LatentRunResults$Intermediary_OLSCoef,LatentRunResults$Intermediary_BootIndex,mean)[-1] )),
       "OLSTstat" = m1_/se1_,
       
       "IVRegCoef" = (m2_ <- tapply(LatentRunResults$Intermediary_IVRegCoef,LatentRunResults$Intermediary_BootIndex,mean)[1]),
       "IVRegSE" = (se2_ <- sd(tapply(LatentRunResults$Intermediary_IVRegCoef,LatentRunResults$Intermediary_BootIndex,mean)[-1] )),
       "IVRegTstat" = m2_/se2_,
       
       "x.est1" = LatentRunResults$Intermediary_x.est1[,1],
       "x.est2" = LatentRunResults$Intermediary_x.est2[,1],
       
       "Corrected_IVRegCoef" = (m3_ <- LatentRunResults$Intermediary_BootInde[1]),
       "Corrected_IVRegSE" = (se3_ <- sd(LatentRunResults$Intermediary_BootIndex[-1])),
       "Corrected_IVRegTstat"  = m3_/se3_,
       "VarEst_split" = VarEst_split,
       "VarEst_split_se" = VarEst_split_se 
  )
  
  return( LatentRunResults ) 
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

