% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpme_DoBootPartition.R
\name{lpme}
\alias{lpme}
\title{lpme}
\usage{
lpme(
  Yobs,
  ObservablesMat,
  ObservablesGroupings = colnames(ObservablesMat),
  MakeObservablesGroupings = F,
  nBoot = 32L,
  nPartition = 10L,
  bootBasis = 1:length(Yobs),
  ReturnIntermediaries = T,
  seed = runif(1, 1, 10000)
)
}
\arguments{
\item{Yobs}{A vector of observed outcome variables}

\item{ObservablesMat}{A matrix of observable indicators used to estimate the latent variable}

\item{ObservablesGroupings}{A vector specifying groupings for the observable indicators. Default is column names of ObservablesMat.}

\item{MakeObservablesGroupings}{Logical. If TRUE, creates dummy variables for each level of the observable indicators. Default is FALSE.}

\item{nBoot}{Integer. Number of bootstrap iterations. Default is 32.}

\item{nPartition}{Integer. Number of partitions for each bootstrap iteration. Default is 10.}

\item{bootBasis}{Vector of indices or grouping variable for stratified bootstrap. Default is 1:length(Yobs).}

\item{ReturnIntermediaries}{Logical. If TRUE, returns intermediate results. Default is TRUE.}

\item{seed}{Random seed for reproducibility. Default is a random integer between 1 and 10000.}
}
\value{
A list containing various estimates and statistics:
\itemize{
  \item OLSCoef: Coefficient from naive OLS regression
  \item OLSSE: Bootstrap standard error of naive OLS coefficient
  \item OLSTstat: T-statistic of naive OLS coefficient
  \item Corrected_OLSCoef: OLS coefficient corrected for measurement error
  \item Corrected_OLSSE: Bootstrap standard error of corrected OLS coefficient
  \item Corrected_OLSTstat: T-statistic of corrected OLS coefficient
  \item Corrected_OLSCoef_alt: Alternative corrected OLS coefficient
  \item Corrected_OLSSE_alt: Bootstrap standard error of alternative corrected OLS coefficient
  \item Corrected_OLSTstat_alt: T-statistic of alternative corrected OLS coefficient
  \item IVRegCoef: Coefficient from instrumental variable regression
  \item IVRegSE: Bootstrap standard error of IV regression coefficient
  \item IVRegTstat: T-statistic of IV regression coefficient
  \item Corrected_IVRegCoef: IV regression coefficient corrected for measurement error
  \item Corrected_IVRegSE: Bootstrap standard error of corrected IV coefficient
  \item Corrected_IVRegTstat: T-statistic of corrected IV coefficient
  \item x.est1: First set of latent variable estimates
  \item x.est2: Second set of latent variable estimates
  \item VarEst_split: Estimated variance of the measurement error
  \item VarEst_split_se: Bootstrap standard error of the estimated variance of the measurement error
}
}
\description{
Implements bootstrapped analysis for latent variable models with measurement error correction
}
\details{
This function implements a bootstrapped latent variable analysis with measurement error correction. 
It performs multiple bootstrap iterations, each with multiple partitions. For each partition, 
it calls the LatentOneRun function to estimate latent variables and apply various correction methods. 
The results are then aggregated across partitions and bootstrap iterations to produce final estimates 
and bootstrap standard errors.
}
\examples{
# Generate some example data
set.seed(123)
Yobs <- rnorm(1000)
ObservablesMat <- as.data.frame( matrix(sample(c(0,1), 1000*10, replace = T), ncol = 10) )

# Run the bootstrapped analysis
results <- lpme(Yobs, ObservablesMat, nBoot = 100, nPartition = 5)

# View the corrected IV coefficient and its standard error
print(c(results$Corrected_IVRegCoef, results$Corrected_IVRegSE))

}
