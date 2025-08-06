#!/bin/bash

# Part 1: install necessary apps
## Install correct R version
R_VERSION_REQUIRED=$(grep '"Version":' renv.lock | head -n 1 | sed 's/[^0-9.]*//g') # Extract required R version from renv.lock

R_VERSION=$(R --version | grep "^R version" | sed -E 's/^R version ([0-9]+\.[0-9]+\.[0-9]+).*/\1/')

if [ "$(printf '%s\n' "$R_VERSION" "$R_VERSION_REQUIRED" | sort -V | head -n1)" = "4.3.2" ]; then
  echo "R version is >= $R_VERSION_REQUIRED. Requirement met."
else
  echo "R version is < $R_VERSION_REQUIRED. Please install (https://cran.rstudio.com/) and switch to R >= $R_VERSION_REQUIRED. We recommend using rig (https://github.com/r-lib/rig) to manage R versions."
  exit 1
fi

## Install gfortran, which is required to to compile FORTRAN source code during R package installation
if ! command -v gfortran >/dev/null 2>&1; then
  echo "GFortran not found. Please install GFortran through https://fortran-lang.org/learn/os_setup/install_gfortran/. For MacOS user, refer to https://mac.r-project.org/tools/."
  exit 1
fi

## Install nc-config (NetCDF configuration script), which is required to install R package 'ncdf4'
if ! command -v nc-config >/dev/null 2>&1; then
  echo "NetCDF (nc-config) not found. Please install NetCDF library through https://github.com/Unidata/netcdf-c?tab=readme-ov-file"
  exit 1
fi

# Part 2: restore R environment
## Install the `renv` package, and restore the project-specific packages specified in `renv.lock`
Rscript setup.R
