# lpmec 1.1.4

## Package Rename
* Renamed package from `lpme` to `lpmec` (Latent Predictor Measurement Error Correction).
* The name "lpme" was already taken by a different archived CRAN package (by Zhou & Huang).
* Function names (`lpme()`, `lpme_onerun()`) remain unchanged for backward compatibility.

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
* Updated version year in CITATION file.

## Previous Changes
* Added CITATION file for proper citation information.
* Updated DESCRIPTION with explicit Author field.
* Minor documentation updates.
