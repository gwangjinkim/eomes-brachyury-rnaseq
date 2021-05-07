
"
binquire <- function(packagename, ...) {
  pckgname <- toString(substitute(packagename))
  eval(substitute(
    if (!require(pckgname)) {
      source('https://bioconductor.org/biocLite.R')
      biocLite(pckgname, ...)
      require(pckgname)
    }
  ))
}
" # old version for documentation


binquire <- function(packagename, ...) {
  pckgname <- toString(substitute(packagename))
  eval(substitute(
    if (!require(pckgname)) {
      BiocManager::install(pckgname, ...)
      require(pckgname)
    }
  ))
}

inquire <- function(packagename, ...) {
  pckgname <- toString(substitute(packagename))
  eval(substitute(
    if (!require(pckgname)) {
      install.packages(pckgname, ...)
      require(pckgname)
    }
  ))
}

