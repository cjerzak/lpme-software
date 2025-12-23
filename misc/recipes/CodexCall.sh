
----------------------------------------------------
# Codex Build
----------------------------------------------------

# 1. Install system-level R, git, and build tools
apt-get update \
  && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
       r-base \
       r-base-dev \
       libcurl4-openssl-dev \
       libssl-dev \
       libxml2-dev \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

#------------------------------------------------------------------------------
# 2. All required CRAN packages via Debian/Ubuntu binary builds
#    (stats ships with base R, so no package needed)
#------------------------------------------------------------------------------
apt-get update \
  && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
       r-cran-remotes \
       r-cran-testthat \
       r-cran-reticulate \
       r-cran-pscl \
       r-cran-aer \
       r-cran-sandwich \
       r-cran-mvtnorm \
       r-cran-amelia \
       r-cran-rcpp \
       r-cran-rcpparmadillo \
       r-cran-gtools \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*


# Quickly install other binaries 
BIN_ROOT="https://raw.githubusercontent.com/cjerzak/lpmec-software/main/misc/binaries"

# Tarballs you produced in the builder container
for f in emIRT_0.0.14_R_x86_64-pc-linux-gnu.tar.gz
do
  curl -sSL -o "/tmp/$f" "$BIN_ROOT/$f"
done

# First library path where R looks for user packages
LIB="$(Rscript -e 'cat(.libPaths()[1])')"
mkdir -p "$LIB"

# Untar each package directly into the library path (no R CMD INSTALL)
for f in /tmp/*.tar.gz; do
  tar -xzf "$f" -C "$LIB"
done
rm /tmp/*.tar.gz

# main install 
Rscript -e "
install.packages('sensemakr');
remotes::install_github(
  'cjerzak/lpmec-software',
  subdir = 'lpmec',
  dependencies = FALSE,
  build_vignettes = FALSE   # speeds up CI builds
)"


