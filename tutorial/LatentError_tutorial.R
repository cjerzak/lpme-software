#!/usr/bin/env Rscript
{
  library(latenterror)
  
  # Generate data 
  Yobs <- rnorm(1000)
  ObservablesMat <- as.data.frame( matrix(sample(c(0,1), 1000*10, replace = T), ncol = 10) )
  
  # One run
  OneRun <- LatentOneRun(Yobs, ObservablesMat, 
               MakeObservablesGroupings = FALSE, seed = runif(1, 1, 10000))
  
  # Latent error correction method, with partitioning and bootstrap 
  FullRun <- LatentRun(Yobs, ObservablesMat, ObservablesGroupings = colnames(ObservablesMat),
            MakeObservablesGroupings = FALSE, nBoot = 32L, nPartition = 10L,
            bootBasis = 1:length(Yobs), ReturnIntermediaries = TRUE,
            seed = runif(1, 1, 10000)) 
}
