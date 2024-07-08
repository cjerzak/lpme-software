# `latenterror`:  An `R` Packge for Measurement Error Corrections Under Identification Restrictions
[**Installation**](#installation)
| [**Key Functions**](#keyfxns)
| [**References**](#references)

`latenterror` is an R package that provides tools for analyzing latent variable models with measurement error correction, using bootstrapping techniques for robust estimation and inference.

# Package Installation<a id="installation"></a>
Within an `R` session, you can install the development version of `latenterror` from GitHub with:
```
install.packages("devtools")
devtools::install_github("cjerzak/latenterror")
```

# Key Functions<a id="keyfxns"></a>
## `LatentOneRun`
`LatentOneRun` performs a single run of latent variable analysis with measurement error correction.
```
LatentOneRun(Yobs, ObservablesMat, ObservablesGroupings = colnames(ObservablesMat),
MakeObservablesGroupings = FALSE, seed = runif(1, 1, 10000))
```

## `LatentRun`
`LatentRun` implements a bootstrapped analysis for latent variable models with measurement error correction.
```
LatentRun(Yobs, ObservablesMat, ObservablesGroupings = colnames(ObservablesMat),
MakeObservablesGroupings = FALSE, nBoot = 32L, nPartition = 10L,
bootBasis = 1:length(Yobs), ReturnIntermediaries = TRUE,
seed = runif(1, 1, 10000))
```

### Example
Here's an example of how to use the main `LatentRun` function:
```
library(latenterror)

# Generate some example data
set.seed(123)
Yobs <- rnorm(1000)
ObservablesMat <- matrix(rnorm(10000), ncol = 10)

# Run the bootstrapped analysis
results <- latenterror::LatentRun(Yobs, ObservablesMat, nBoot = 100, nPartition = 5)

#View the corrected IV coefficient and its standard error
print(c(results$Corrected_IVRegCoef, results$Corrected_IVRegSE))
```

# Contributing
Contributions to latenterror are welcome! Feel free to submit a [pull request](https://github.com/cjerzak/latenterror-software/pulls) or open an [issue](https://github.com/cjerzak/latenterror/issues).

# Acknowledgements 
We thank [Jeff Lewis](https://polisci.ucla.edu/person/jeffrey-b-lewis/), [Umberto Mignozzetti](https://umbertomig.com/), [Aaron Pancost](https://sites.google.com/site/aaronpancost/), [Erik Snowberg](https://eriksnowberg.com/), [Chris Tausanovitch](https://ctausanovitch.com/), and participants of a panel at an MPSA panel for very helpful comments. We thank [Major Valls](https://www.linkedin.com/in/major-valls-39b6b9229/) for excellent research assistance.

# References<a id="references"></a>
[Connor T. Jerzak](https://github.com/cjerzak), [Stephen A. Jessee](https://github.com/sjessee). Measurement Error in Latent Predictors: The Role of Identification Restrictions *Working paper to be released summer 2024!*
