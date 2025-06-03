

----------------------------------------------------
# Docker Build (place x86 tarballs to binaries folder in GitHub)
----------------------------------------------------

docker run --platform=linux/amd64 --rm \
  -v "$(pwd)/binaries:/binaries" \
  rocker/r-ver:4.4.0 bash -exc "
    set -euo pipefail

    ## 1) Install system development libraries (no R packages here)
    echo '🔧 Installing system development libraries…'
    apt-get update -qq
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      build-essential \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev

    ## 2) Make sure R can install pscl (and any other pure‐R deps) from CRAN
    echo '📦 Installing pscl from CRAN…'
    Rscript -e 'install.packages(
      \"pscl\",
      repos = \"https://cloud.r-project.org\",
      dependencies = TRUE
    )'

    ## 3) Switch into /binaries so download.packages writes there
    cd /binaries

    ################################################################################
    ## 4) Build Rcpp → “Rcpp_<ver>_R_x86_64-pc-linux-gnu.tar.gz”
    echo '📦 Downloading Rcpp source…'
    Rscript -e 'download.packages(
      \"Rcpp\",
      destdir      = \"/binaries\",
      type         = \"source\",
      repos        = \"https://cloud.r-project.org\"
    )'
    echo '⚙️  Building Rcpp binary…'
    R CMD INSTALL --build /binaries/Rcpp_*.tar.gz

    ################################################################################
    ## 5) Build RcppArmadillo → “RcppArmadillo_<ver>_R_x86_64-pc-linux-gnu.tar.gz”
    echo '📦 Downloading RcppArmadillo source…'
    Rscript -e 'download.packages(
      \"RcppArmadillo\",
      destdir      = \"/binaries\",
      type         = \"source\",
      repos        = \"https://cloud.r-project.org\"
    )'
    echo '⚙️  Building RcppArmadillo binary…'
    R CMD INSTALL --build /binaries/RcppArmadillo_*.tar.gz

    ################################################################################
    ## 6) Build emIRT → “emIRT_<ver>_R_x86_64-pc-linux-gnu.tar.gz”
    echo '📦 Downloading emIRT source…'
    Rscript -e 'download.packages(
      \"emIRT\",
      destdir      = \"/binaries\",
      type         = \"source\",
      repos        = \"https://cloud.r-project.org\"
    )'
    echo '⚙️  Building emIRT binary…'
    R CMD INSTALL --build /binaries/emIRT_*.tar.gz

    ################################################################################
    ## 7) List out exactly what landed in /binaries
    echo '✅ Built binary tarballs:'
    ls -l /binaries
  "
