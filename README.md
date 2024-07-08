# `latenterror`:  An `R` Packge for Measurement Error Corrections under Identification Restriction

`latenterror` is an R package that provides tools for analyzing latent variable models with measurement error correction, using bootstrapping techniques for robust estimation and inference.

#Installation
You can install the development version of latenterror from GitHub with:
```
install.packages("devtools")
devtools::install_github("cjerzak/latenterror")
```

Key Functions
`LatentOneRun`
This function performs a single run of latent variable analysis with measurement error correction.
```
LatentOneRun(Yobs, ObservablesMat, ObservablesGroupings = colnames(ObservablesMat),
MakeObservablesGroupings = FALSE, seed = runif(1, 1, 10000))
```
LatentRun
This function implements a bootstrapped analysis for latent variable models with measurement error correction.
```
LatentRun(Yobs, ObservablesMat, ObservablesGroupings = colnames(ObservablesMat),
MakeObservablesGroupings = FALSE, nBoot = 32L, nPartition = 10L,
bootBasis = 1:length(Yobs), ReturnIntermediaries = TRUE,
seed = runif(1, 1, 10000))
```
Example
Here's an example of how to use the LatentRun function:
```
library(latenterror)

# Generate some example data
set.seed(123)
Yobs <- rnorm(1000)
ObservablesMat <- matrix(rnorm(10000), ncol = 10)

# Run the bootstrapped analysis
results <- latenterror::LatentRun(Yobs, ObservablesMat, nBoot = 100, nPartition = 5)
```
View the corrected IV coefficient and its standard error
```
print(c(results$Corrected_IVRegCoef, results$Corrected_IVRegSE))
```

# Contributing
Contributions to latenterror are welcome. Please feel free to submit a Pull Request or open an (issue)[https://github.com/cjerzak/latenterror/issues]. 
