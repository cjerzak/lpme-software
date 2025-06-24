# `lpme`: R Package for Dealing with Latent Predictor Measurement Error Under Identification Restrictions
[**Installation**](#installation)
| [**Key Functions**](#keyfxns)
| [**References**](#references)

[<img src="https://img.shields.io/badge/Demo-View%20Demo-blue" alt="Demo Button">](https://github.com/cjerzak/lpme-software/blob/main/lpme/vignettes/IntroVignette.Rmd)
[![CI](https://github.com/cjerzak/lpme-software/actions/workflows/ci.yml/badge.svg)](https://github.com/cjerzak/lpme-software/actions/workflows/ci.yml)

`lpme` is an R package that provides tools for analyzing latent variable models with measurement error correction, using bootstrapping techniques for inference.

# Package Installation<a id="installation"></a>
Within an `R` session, you can install the development version of `lpme` from GitHub with:
```
# Install from GitHub
# install.packages("devtools") 
devtools::install_github("cjerzak/lpme-software/lpme")
```

# Key Functions<a id="keyfxns"></a>
## `lpme_onerun`
`lpme_onerun` performs a single run of latent variable analysis with measurement error correction (no bootstrapping; 1 split sample partition): 
```
# Generate data 
Yobs <- rnorm(1000)
ObservablesMat <- matrix(sample(c(0,1), 1000*10, replace = TRUE), ncol = 10)

# One run of latent error correction method 
lpme::lpme_onerun(Y = Yobs, 
                  observables = ObservablesMat)
```

## `lpme`
`lpme` implements a bootstrap analysis for latent variable models with measurement error correction. We average over `n_partition` split sample partitions. 
```
# Generate data
Yobs <- rnorm(1000)
ObservablesMat <- matrix(sample(c(0,1), 1000*10, replace = TRUE), ncol = 10)

# Latent error correction method, with partitioning and bootstrap
results <- lpme::lpme(
    Y = Yobs,
    observables = ObservablesMat,
    n_boot = 32L,
    n_partition = 10L
)

# View the corrected IV coefficient and its standard error
print(results)
```

# Contributing
Contributions to lpme are welcome! Feel free to submit a [pull request](https://github.com/cjerzak/lpme-software/pulls) or open an [issue](https://github.com/cjerzak/lpme-software/issues).

# Acknowledgements 
We thank [Guilherme Duarte](https://duarteguilherme.github.io/), [Jeff Lewis](https://polisci.ucla.edu/person/jeffrey-b-lewis/), [Umberto Mignozzetti](https://umbertomig.com/), [Aaron Pancost](https://sites.google.com/site/aaronpancost/), [Erik Snowberg](https://eriksnowberg.com/), [Chris Tausanovitch](https://ctausanovitch.com/), and participants of an MPSA panel for very helpful comments. We thank [Major Valls](https://www.linkedin.com/in/major-valls-39b6b9229/) for excellent research assistance.

# References<a id="references"></a>
[Connor T. Jerzak](https://github.com/cjerzak), [Stephen A. Jessee](https://github.com/sjessee). Measurement Error in Latent Predictors: The Role of Identification Restrictions. _Working Paper._

[<img src="https://connorjerzak.com/wp-content/uploads/2024/08/LatentErrorViz.png">](https://connorjerzak.com/)
