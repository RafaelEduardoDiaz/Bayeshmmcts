#' Bayeshmmcts: Fit hidden Markov Models via RStan.
#'
#' @description Fast and flexible tool for Bayesian estimation of hidden Markov models using the 'rstan' package,
#'              which provides the R interface to the Stan C++ library. At the moment only the PHMM and ZIPHMM are implemented.
#' @docType package
#' @name Bayeshmmcts-package
#' @aliases Bayeshmmcts
#' @useDynLib Bayeshmmcts, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'
NULL
