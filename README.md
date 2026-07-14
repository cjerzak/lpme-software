# `lpmec`: Measurement Error Analysis and Correction under Identification Restrictions 
[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [<img src="https://img.shields.io/badge/Demo-View%20Demo-blue" alt="Demo Button">](https://github.com/cjerzak/lpmec-software/blob/main/lpmec/vignettes/IntroVignette.Rmd) [<img src="https://img.shields.io/badge/CRAN-View%20on%20CRAN-green" alt="CRAN Button">](https://cran.r-project.org/package=lpmec)

[**Installation**](#installation) | [**Key Functions**](#keyfxns) | [**References**](#references)

`lpmec` is an R package that provides tools for analyzing latent variable models with measurement error correction, using bootstrapping techniques for inference.

## Overview
Measurement error in latent predictors (e.g., ideology scores from survey responses, ability measures from test items) causes attenuation bias in regression coefficients—systematically biasing estimates toward zero.

`lpmec` implements split-sample instrumental variables and OLS correction methods to obtain consistent estimates when your predictor is measured with error.

**Key features:**
- Split-sample IV estimation for latent predictors
- Multiple estimation methods: EM (emIRT), PCA, MCMC (pscl/NumPyro)
- Bootstrap inference with confidence intervals
- Sensitivity analysis via sensemakr integration

# Package Installation<a id="installation"></a>
Within an `R` session, you can install the development version of `lpmec` from GitHub with:
```
# Install from GitHub
# install.packages("devtools")
devtools::install_github("cjerzak/lpmec-software", subdir = "lpmec")
```

## Quick Start
```r
library(lpmec)
data(KnowledgeVoteDuty)

# Run correction with bootstrap
results <- lpmec(Y = KnowledgeVoteDuty$voteduty,
                 observables = as.matrix(KnowledgeVoteDuty[, 2:5]),
                 n_boot = 50)

print(results)
# Latent Predictor Measurement Error Correction (LPMEC) Model Results
# -------------------------------------------------------------------
# Uncorrected Coefficient (OLS): X.XXX (SE: X.XXX)
# Corrected Coefficient: X.XXX (SE: X.XXX)
```

## Results

The correlation-corrected estimates computed here by `lpmec` closely track the performance of the computationally intensive full Bayesian joint estimation, providing good results regardless of the sample size (N) or number of indicators (M). 

### Bias \& Variance 

<a href="https://arxiv.org/abs/2507.22218/#gh-light-mode-only">
  <img src="https://connorjerzak.com/wp-content/uploads/2026/01/p2.webp#gh-light-mode-only" alt="Figure – light" width="500">
</a>

<a href="https://arxiv.org/abs/2507.22218/#gh-dark-mode-only">
  <img src="https://connorjerzak.com/wp-content/uploads/2026/01/p2_dark.webp#gh-dark-mode-only" alt="Figure – dark" width="500">
</a>

### Coverage 

<a href="https://arxiv.org/abs/2507.22218/#gh-light-mode-only">
  <img src="https://connorjerzak.com/wp-content/uploads/2026/01/p5.webp#gh-light-mode-only" alt="Figure – light" width="350">
</a>

<a href="https://arxiv.org/abs/2507.22218/#gh-dark-mode-only">
  <img src="https://connorjerzak.com/wp-content/uploads/2026/01/p5_dark.webp#gh-dark-mode-only" alt="Figure – dark" width="350">
</a>




## Optional: NumPyro Backend for MCMC

For advanced MCMC estimation, `lpmec` supports the NumPyro backend via Python. NumPyro leverages JAX’s automatic differentiation and JIT/XLA compilation (including optional GPU/TPU execution), which can make HMC/NUTS sampling faster and easier to parallelize than CPU-only MCMC backends. This is optional; the default `pscl` backend works without any Python setup.

To use the NumPyro backend:

1. **Install Python dependencies** using the built-in setup function:
```r
library(lpmec)
build_backend()  # Creates conda environment with JAX and NumPyro
```

2. **Alternative: Manual installation** if you prefer to manage your own Python environment:
```bash
# Create and activate a conda environment
conda create -n lpmec python=3.10
conda activate lpmec

# Install JAX and NumPyro
pip install jax jaxlib numpyro
```

3. **Use the NumPyro backend** in your analysis:
```r
results <- lpmec_onerun(
  Y = Yobs,
  observables = ObservablesMat,
  estimation_method = "mcmc",
  mcmc_control = list(backend = "numpyro")
)
```

Note: The NumPyro backend requires a working Python installation accessible via the `reticulate` package.

# Key Functions<a id="keyfxns"></a>
## `lpmec_onerun`
`lpmec_onerun` performs a single run of latent variable analysis with measurement error correction (no bootstrapping; 1 split sample partition):
```
# Generate data
Yobs <- rnorm(1000)
ObservablesMat <- matrix(sample(c(0,1), 1000*10, replace = TRUE), ncol = 10)

# One run of latent error correction method
lpmec::lpmec_onerun(Y = Yobs,
                    observables = ObservablesMat)
```

## `lpmec`
`lpmec` implements latent variable models with measurement error correction. It fits the original sample, optionally adds bootstrap resampling via `n_boot`, and averages over `n_partition` split sample partitions. Set `n_boot = 0` to disable bootstrap resampling while still aggregating original-sample partitions.
```
# Generate data
Yobs <- rnorm(1000)
ObservablesMat <- matrix(sample(c(0,1), 1000*10, replace = TRUE), ncol = 10)

# Latent error correction method, with partitioning and bootstrap
results <- lpmec::lpmec(
    Y = Yobs,
    observables = ObservablesMat,
    n_boot = 32L,
    n_partition = 10L
)

# View the corrected IV coefficient and its standard error
print(results)
```

`lpmec` ships with a small example dataset, `KnowledgeVoteDuty`, drawn from the 2024 American National Election Study. Load it with `data(KnowledgeVoteDuty)` to try the package on real survey responses.

## Estimation Methods
| Method | Description | Backend |
|--------|-------------|---------|
| `em` (default) | EM algorithm | emIRT |
| `pca` | Principal component | base R |
| `averaging` | Simple row means | base R |
| `mcmc` | Full Bayesian MCMC | pscl or NumPyro |
| `mcmc_joint` | Joint Bayesian latent-variable and outcome model | NumPyro |
| `mcmc_joint2` | NumPyro mixed factor-analysis benchmark | NumPyro |
| `mcmc_overimputation` | MCMC with over-imputation correction | pscl or NumPyro |
| `custom` | User-provided function | - |

## Documentation
- [Introduction Vignette](https://github.com/cjerzak/lpmec-software/blob/main/lpmec/vignettes/IntroVignette.Rmd) - Simulated-data walkthrough of the basic workflow
- Run `?lpmec` after installation for function documentation

# Contributing
Contributions to lpmec are welcome! Feel free to submit a [pull request](https://github.com/cjerzak/lpmec-software/pulls) or open an [issue](https://github.com/cjerzak/lpmec-software/issues).

# Acknowledgements 
We thank [Guilherme Duarte](https://duarteguilherme.github.io/), [Jeff Lewis](https://polisci.ucla.edu/person/jeffrey-b-lewis/), [Umberto Mignozzetti](https://umbertomig.com/), [Aaron Pancost](https://sites.google.com/site/aaronpancost/), [Erik Snowberg](https://eriksnowberg.com/), [Chris Tausanovitch](https://ctausanovitch.com/), and participants of an MPSA panel for very helpful comments. We thank [Major Valls](https://www.linkedin.com/in/major-valls-39b6b9229/) for excellent research assistance.

# Citation<a id="references"></a>
[Connor T. Jerzak](https://github.com/cjerzak), [Stephen
A. Jessee](https://github.com/sjessee). Correcting Bias When Using Latent Regressors. Forthcoming at *Political Analysis* (pending final replication), 2026. [[PDF]](https://arxiv.org/pdf/2507.22218)

```bibtex
@article{jerzak2026attenuationbiaslatent,
      title={Correcting Bias When Using Latent Regressors}, 
      author={Jerzak, Connor T. and Jessee, Stephen A.},
      journal={Political Analysis},
      year={2026},
      note={Forthcoming pending final replication},
      url={https://arxiv.org/abs/2507.22218},
      keywords={latent variables, measurement error, attenuation bias, identification restrictions, method of composition}
}
```

<!---
[<img src="https://connorjerzak.com/wp-content/uploads/2024/08/LatentErrorViz.png">](https://connorjerzak.com/)
-->


<!---
# Happy lpmec-ing! 

```
                        True Values
                     [Mean 0, Var = 1]

                    o   oo ooooo oo   o
                   /   /  /  |  \  \   \
                  /   /  /   |   \  \   \
                 /   /  /    |    \  \   \
                /   /  /     |     \  \   \
               /   /  /      |      \  \   \
              v   v  v       v       v  v   v

             o    o  o     ooooo     o  o    o
             
                      With Noise
               [Mean 0, Var = 1 + σ_u²]

             o    o  o     ooooo     o  o    o
              \    \  \      |      /  /    /
               \    \  \     |     /  /    /
                \    \  \    |    /  /    /
                 \    \  \   |   /  /    /
                  \    \  \  |  /  /    /
                   v    v  v v v  v    v
                     ●   ●●●●●●●   ●

        With Noise & Data-dependent Normalization
                    [Mean 0, Var = 1]


       |-------|-------|-------|-------|-------|-------|
      -3      -2      -1       0       1       2       3
```
-->
