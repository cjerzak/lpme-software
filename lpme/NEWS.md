# lpme 0.1.1

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
