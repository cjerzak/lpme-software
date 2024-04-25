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

LatentErrorRun <- function(Yobs,
                           AnalyzeFactors,
                           factorizePolicies = F){
  items.split1 <- sample(1:length(AnalyzeFactors), size = floor(length(AnalyzeFactors)/2), replace=F)
  items.split2 <- (1:length(AnalyzeFactors))[! (1:length(AnalyzeFactors) %in% items.split1)]
  for(split_ in c("", 1, 2)){
    if(split_ == ""){ items.split_ <- AnalyzeFactors[1:length(AnalyzeFactors)] }
    if(split_ == "1"){ items.split_ <- AnalyzeFactors[items.split1] }
    if(split_ == "2"){ items.split_ <- AnalyzeFactors[items.split2] }
    
    # estimating ideal points
    if(factorizePolicies == F){
      POLICIES <- POLICIES[,items.split_]
    }
    if(factorizePolicies == T){
      POLICIES <- do.call(cbind,unlist(apply(POLICIES[,items.split_],2,function(zer){
        list( model.matrix(~0+as.factor(zer))[,-1] )}),recursive = F))
    }
    x_init <- apply( POLICIES[,items.split_], 1, function(x){ mean(f2n(x), na.rm=T)})
    rc_ <- convertRC( RC.SIM_ <- rollcall(POLICIES) )
    p_ <- makePriors(rc_$n, rc_$m, 1)
    s_ <- list("alpha" = matrix(rnorm(ncol(POLICIES))),
               "beta" = matrix(rnorm(ncol(POLICIES))),
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
  
  list("SimpleRegCoef" = coef(simpleReg)[2],
       "SimpleRegTstat" = coef(summary(simpleReg))[2,3],
       "IVRegTstat" = coef(summary(IVStage2))[2,3],
       "IVRegCoef" = coef(IVStage2)[2], 
       "Corrected_IVRegCoef" = coef(IVStage2)[2] / sqrt(1 + var(my_df_$x.est1 - my_df_$x.est2)/2),
       "Corrected_IVRegTstat" =  Corrected_IVRegCoef[monti_i] / coef(summary(IVStage2))[2,2],
       "VarEst_split" = var(my_df_$x.est1 - my_df_$x.est2) / 2)
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

