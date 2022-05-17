# python 'sympy' module I want to use in my package
sympy <- NULL

.onLoad <- function(libname, pkgname) {
  # delay load foo module (will only be loaded when accessed via $)
  scipy <<- reticulate::import("sympy", delay_load = TRUE)
}
