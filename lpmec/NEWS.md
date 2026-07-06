# lpmec 1.2.0

## Panel and Fixed-Effects Designs
* Added `lpmec_panel_onerun()` and `lpmec_panel()` for measurement error
  correction under panel designs (`"within"`, `"twoway"`, `"difference"`,
  `"pooled"`): NA-safe design transforms (iterated demeaning for unbalanced
  two-way panels; exact-gap k-period differencing), cluster-robust naive OLS,
  within- and cross-measure split-IV, and the design-aware correction algebra
  (corrected OLS is `b * sqrt(rho_pooled) / rho_design`; corrected split-IV
  multiplies by `sqrt(rho_pooled)` -- IV multiplies, OLS divides).
* `lpmec_panel()` reports cluster (unit) bootstrap uncertainty with fresh
  pseudo-ids per draw, re-running the whole pipeline (scoring, sign
  alignment, transforms, reliabilities, regressions, corrections) on every
  (bootstrap, partition) pair.

## Reliability Bounds
* Added `lpmec_reliability_bounds()`: per (measure, design) split
  correlations and cross-measure triad reliabilities
  (`rho*_m = r_ml * r_mk / r_lk`, requiring 3+ measures), with the implied
  identification interval `[rho_lo, rho_hi]` and an optional cluster
  bootstrap. Reliabilities below `min_reliability` and correlations on fewer
  than `min_cor_n` complete pairs are reported as `NA` rather than divided
  out.

## Latent Moderators
* Added `lpmec_moderator_onerun()` and `lpmec_moderator()` for
  treatment-by-latent-moderator interactions in experiments: the naive
  interaction coefficient is divided by the square root of the
  Spearman-Brown as-used score reliability
  (`rho_score = M * rbar / (1 + (M - 1) * rbar)`), with per-measure triad
  corrections when 3+ measures are supplied and a row (or stratified)
  bootstrap that propagates reliability-estimation uncertainty into the
  corrected-interaction intervals.

## Measure Inputs
* The new functions share a tri-source measure-input model merged by name:
  `observables` (raw item batteries scored via `lpmec_onerun()`),
  `split_scores` (n x 2 half-score matrices, required for measures with
  fewer than 4 items), and `scores` (pre-computed full scores). All scores
  are pooled z-scored and sign-aligned before any correlations are taken.

## S3 Methods
* Added `print()`/`summary()` methods for `lpmec_panel_onerun`,
  `lpmec_panel`, `lpmec_moderator_onerun`, and `lpmec_moderator`,
  `print.lpmec_reliability_bounds()`, and `plot()` methods for `lpmec_panel`
  and `lpmec_moderator`.

## Documentation
* Vignette gains three sections: panel/fixed-effects corrections,
  reliability bounds with multiple measures, and latent moderators in
  experiments.

## Bug Fixes
* Fixed `.lpmec_aggregate_by_boot()` dropping the matrix structure for
  single-column inputs, which broke bootstrap aggregation for
  `lpmec_multivariate()` runs with a single latent predictor.

# lpmec 1.1.4

## Bootstrap Inference
* Added m-out-of-n and subsampling uncertainty for nonsmooth finite-partition
  median aggregation. The implementation fixes the realized partition set
  across resamples, reruns the full latent measurement and correction pipeline
  within each resample, and reports root-scaled standard errors and confidence
  intervals for m < n designs.

## Package Rename
* Renamed package from `lpme` to `lpmec` (Latent Predictor Measurement Error Correction).
* The name "lpme" was already taken by a different archived CRAN package (by Zhou & Huang).
* Main functions are `lpmec()` and `lpmec_onerun()`.

## CRAN Resubmission
* Removed dependency on archived package `decon`.
* Version bump from 0.1.1 to 1.1.4 (continuing from archived CRAN version 1.1.3).

## CRAN Preparation
* Fixed vignette to use correct `estimation_method` values ("em" instead of "emIRT").
* Added `.Rbuildignore` to exclude development files from package build.
* Replaced `eval(parse())` pattern with direct assignment for CRAN compliance.
* Changed `F`/`T` to `FALSE`/`TRUE` throughout codebase.
* Wrapped examples in `\donttest{}` to avoid CRAN timeout issues.
* Changed `conda_env_required` default to `FALSE` for CRAN compatibility.
* Updated CITATION file to use modern `bibentry()` format.
* Added `skip_on_cran()` guards to test files.
* Added internal function documentation with `@noRd` tags.

## Documentation
* Fixed documentation to correctly indicate `pscl` as the default MCMC backend.
* Refreshed public documentation for advanced estimation methods, MCMC controls,
  return-value fields, and current arXiv links.
* Updated version year in CITATION file.

## Previous Changes
* Added CITATION file for proper citation information.
* Updated DESCRIPTION with explicit Author field.
* Minor documentation updates.
