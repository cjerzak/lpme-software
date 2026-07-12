# lpmec 1.3.0

## Correction Fixes (results-changing)
* **Split-based corrected OLS now uses score-scale (Spearman-Brown)
  reliabilities.** In `split_scores` and `observables` modes the regression
  uses the full score (the average of the two z-scored halves, or the
  full-battery estimate) while the raw half-split correlation estimates the
  reliability of a *half* score, so the previous correction
  `b * sqrt(r_pooled) / r_design` overcorrected (by tens of percent under
  aggressive design transforms). The split variant now applies
  `b * sqrt(SB(r_pooled)) / SB(r_design)` with `SB(r) = 2r / (1 + r)`,
  matching Assumption 3 of the accompanying paper and agreeing
  asymptotically with the corrected split-IV. New fields
  `split_rho_score` / `design_split_rho_score` report the stepped-up
  values; the raw correlations are still returned. The triad, pair, and
  split-IV corrections were already correctly scaled and are unchanged, as
  is the two-measure `scores`-mode path (the V2 reference pipeline).
* **Cross-measure IV corrections are now target-oriented (Proposition
  3c).** With 3+ measures, each `m_by_l` coefficient is multiplied by
  `sqrt` of target `m`'s pooled triad reliability instead of `sqrt` of the
  pairwise correlation `r_ml` (consistent only under equal reliabilities).
  The two-measure parallel-pair fallback is unchanged.

## New Methods
* **Latent outcomes (Proposition 4).** `lpmec_panel_onerun()` and
  `lpmec_panel()` accept `Y_split_scores` (an n x 2 matrix of outcome half
  scores): the outcome reliability `rho_Y` is estimated by Spearman-Brown
  step-up and every corrected estimator is additionally divided by
  `sqrt(rho_Y)`. `Y` becomes optional (the outcome score is then built from
  the halves), and `unit` may be `NULL` for `design = "pooled"`, covering
  the pure cross-sectional case.
* **Design-local scale (Proposition 2a).** The panel functions report
  `sd_design_x`, the naive local slope `ols_coef_local` (effect per SD of
  the design-transformed latent trait), and corrected local coefficients
  `corrected_ols_coef_local(_split/_triad/_pair)` with their sensitivity
  range, satisfying the exact local-pooled bridge identity.

## Compatibility
* `Y_split_scores` is inserted after `covariates` in the signatures of
  `lpmec_panel_onerun()` and `lpmec_panel()`; calls that passed `design`
  or later arguments positionally must switch to named arguments.

## Reliability Diagnostics
* `lpmec_reliability_bounds()` and the panel reliability tables gain a
  `rho_split` column (Spearman-Brown step-up of the split correlation to
  the full-score scale); `rho_lo` / `rho_hi` now span `{triad, rho_split}`
  so both candidates estimate the same full-score target.
* Documentation reframed per the paper's Proposition 7: the
  `[rho_lo, rho_hi]` interval is a *sensitivity range*, not a
  partial-identification interval, unless the corresponding orthogonality
  or ratio condition is maintained; directional claims about the split and
  triad candidates are now stated with their conditions.
* Proposition numbering in the documentation updated to the current
  manuscript (design correction = Prop 3, latent controls = Prop 5,
  latent moderators = Prop 6, triangulation = Prop 7).

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
