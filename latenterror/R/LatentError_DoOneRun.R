#!/usr/bin/env Rscript
#' MainFunction
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
#' #x <- MainFunction(x)
#' @export

#' @import Rfast
#' @md

LatentOneRun <- function(Yobs,
                         ObservablesMat, 
                         ObservablesGroupings = colnames(ObservablesMat),
                         MakeObservablesGroupings = F){
  items.split1_names <- sample(unique(ObservablesGroupings), size = floor(length(unique(ObservablesGroupings))/2), replace=F)
  items.split2_names <- unique(ObservablesGroupings)[! (ObservablesGroupings %in% items.split1_names)]
  for(split_ in c("", 1, 2)){
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
    x_init <- apply( ObservablesMat_, 1, function(x){ mean(f2n(x), na.rm=T)})
    rc_ <- convertRC( RC.SIM_ <- rollcall(ObservablesMat_) )
    p_ <- makePriors(rc_$n, rc_$m, 1)
    s_ <- list("alpha" = matrix(rnorm(ncol(ObservablesMat_))),
               "beta" = matrix(rnorm(ncol(ObservablesMat_))),
               "x" = matrix(x_init))
    
    # fixing in case directions of s1 and s2 x starts are flipped
    if(split_ %in% c("1","2")){ if(cor(s_$x, s_past$x) < 0){ s_$x <- -s_$x; s_$beta <- -s_$beta } }
    lout.sim_ <- binIRT(.rc = rc_,
                        .starts = s_, .priors = p_, .control={
                          list(threads=1, verbose=FALSE, thresh=1e-6) },
                        .anchor_subject = which.max(x_init)) # set direction
    x.est_ <- scale(lout.sim_$means$x); s_past <- s_
    if(cor(x.est_, Yobs, use="p") < 0){ x.est_ <- -x.est_ }
    eval(parse(text = sprintf("x.est%s <- x.est_", split_)))
  }
  
  # simple linear reg 
  simpleReg <- lm(Yobs ~ x.est)
  
  # stage 1 results
  IVStage1 <- lm(x.est2 ~ x.est1)
  
  # save baseline IV results
  IVStage2 <- AER::ivreg(Yobs ~ x.est2 | x.est1)
  
  list("OLSCoef" = coef(summary(simpleReg))[1,2],
       "OLSSE" = coef(summary(simpleReg))[2,2],
       "OLSTstat" = coef(summary(simpleReg))[2,3],
       
       "IVRegCoef" = coef(summary(IVStage2))[2,1], 
       "IVRegSE" = coef(summary(IVStage2))[2,2],
       "IVRegTstat" = coef(summary(IVStage2))[2,3],
       
       "x.est1" = x.est1,
       "x.est2" = x.est2,
       
       "Corrected_IVRegCoef" = (Corrected_IVRegCoef <- (coef(IVStage2)[2] / sqrt(1 + var(x.est1 - x.est2)/2))),
       "Corrected_IVRegTstat" =  Corrected_IVRegCoef / coef(summary(IVStage2))[2,2],
       "VarEst_split" = var(x.est1 - x.est2) / 2 
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

