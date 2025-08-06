source("renv/activate.R") # Activate project, which set renv/library as .libPath
# The .Rprofile file is sourced automatically by R when a new R session starts in the project directory, which automatically source `renv/activate.R`.

code <- readLines("renv/activate.R")
version_line <- grep('version\\s*<-\\s*\\"[0-9.]+\\"', code, value = TRUE)
version <- sub('.*version\\s*<-\\s*\\"([0-9.]+)\\".*', '\\1', version_line)

if (!requireNamespace("renv", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes") # Ensure 'remotes' is available
  }
  remotes::install_version("renv", version = version) # Ensure 'renv' is available
}

rm(code, version_line, version)

# renv::activate()
renv::restore()  # Install packages as specified in renv.lock
