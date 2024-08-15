# `lpme`:  R Packge for Dealing with Latent Predictor Measurement Error Under Identification Restrictions
[**Installation**](#installation)
| [**Key Functions**](#keyfxns)
| [**References**](#references)

`lpme` is an R package that provides tools for analyzing latent variable models with measurement error correction, using bootstrapping techniques for inference.

# Package Installation<a id="installation"></a>
Within an `R` session, you can install the development version of `lpme` from GitHub with:
```
# install.packages("devtools") 
devtools::install_github("cjerzak/lpme/lpme-software")
```

# Key Functions<a id="keyfxns"></a>
## `lpme_OneRun`
`lpme_OneRun` performs a single run of latent variable analysis with measurement error correction (no bootstrapping; 1 split sample partition): 
```
# Generate data 
Yobs <- rnorm(1000)
ObservablesMat <- matrix(sample(c(0,1), 1000*10, replace = T), ncol = 10)

# One run of latent error correction method 
lpme::lpme_OneRun(Yobs, ObservablesMat, 
                  MakeObservablesGroupings = FALSE, seed = runif(1, 1, 10000))
```

## `lpme`
`lpme` implements a bootstrap analysis for latent variable models with measurement error correction. We average over `nPartition` split sample partitions. 
```
# Generate data 
Yobs <- rnorm(1000)
ObservablesMat <- matrix(sample(c(0,1), 1000*10, replace = T), ncol = 10)

# Latent error correction method, with partitioning and bootstrap 
results  <- lpme::lpme(Yobs, ObservablesMat, 
	     MakeObservablesGroupings = FALSE, nBoot = 32L, nPartition = 10L,
	     bootBasis = 1:length(Yobs), ReturnIntermediaries = TRUE,
	     seed = runif(1, 1, 10000)) 
		 
#View the corrected IV coefficient and its standard error
print(c(results$Corrected_IVRegCoef, results$Corrected_IVRegSE))
```

# Contributing
Contributions to latenterror are welcome! Feel free to submit a [pull request](https://github.com/cjerzak/latenterror-software/pulls) or open an [issue](https://github.com/cjerzak/latenterror/issues).

# Acknowledgements 
We thank [Guilherme Duarte](https://duarteguilherme.github.io/), [Jeff Lewis](https://polisci.ucla.edu/person/jeffrey-b-lewis/), [Umberto Mignozzetti](https://umbertomig.com/), [Aaron Pancost](https://sites.google.com/site/aaronpancost/), [Erik Snowberg](https://eriksnowberg.com/), [Chris Tausanovitch](https://ctausanovitch.com/), and participants of a panel at an MPSA panel for very helpful comments. We thank [Major Valls](https://www.linkedin.com/in/major-valls-39b6b9229/) for excellent research assistance.

# References<a id="references"></a>
[Connor T. Jerzak](https://github.com/cjerzak), [Stephen A. Jessee](https://github.com/sjessee). Measurement Error in Latent Predictors: The Role of Identification Restrictions *Working paper to be released summer 2024!*

[<img src="https://connorjerzak.com/wp-content/uploads/2024/08/LatentErrorViz.png" width="420" height="400">](https://connorjerzak.com/)
