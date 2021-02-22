  <!-- badges: start -->
  [![Travis build status](https://travis-ci.com/RafaelEduardoDiaz/Bayeshmmcts.svg?branch=master)](https://travis-ci.com/RafaelEduardoDiaz/Bayeshmmcts)
  <!-- badges: end -->
  
# Bayeshmmcts <a href='https://github.com/RafaelEduardoDiaz/Bayeshmmcts/blob/master/man/figure/hmm.png'><img src='man/figure/hmm.png' align="right" height="139" /></a>

## Overview

Fast and flexible tool for Bayesian estimation of hidden Markov models using the 'rstan' package, 
which provides the R interface to the Stan C++ library. At the moment only the PHMM and ZIPHMM are implemented.

## Verify

``` r
pkgbuild::has_build_tools(debug = TRUE)
Found in Rtools 4.0 installation folder
[1] TRUE
```

check rstan package works

``` r
library(rstan)
example(stan_model, package = "rstan", run.dontrun = TRUE)
```

If not, go to the following link [RStan Getting Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) to prepare the installation.

## Installation

First you need to have installed rstan and devtools. Recommended but not necessary install [rtools](https://cran.r-project.org/bin/windows/Rtools/) according your R version. 
Then install `Bayeshmmcts` development version that is in this github repository.

``` r
# install.packages(c("rstan","devtools","rstantools"))
packageVersion("rstantools") #‘2.1.1’
devtools::install_github("RafaelEduardoDiaz/Bayeshmmcts")
```

## Usage

``` r
library(Bayeshmmcts) 
data("homicides") # load data
Homicides <- ts(data = round(homicides$Rate), start = 1960) # convert to time series

# fit stationary Poisson hidden Markov Model 2 states (clasical theory)
PHMM_2states <- pois.HMM.mle(o = Homicides, m = 2,
                             lambda0 = c(30, 63),
                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE), 
                             stationary = TRUE)
print(PHMM_2states)


$m
[1] 2

$lambda
[1] 29.71593 62.81262

$A
           [,1]       [,2]
[1,] 0.98024196 0.01975804
[2,] 0.06396795 0.93603205

$pi
[1] 0.7640155 0.2359845

$code
[1] 1

$mllk
[1] -201.3239

$AIC
[1] 410.6477

$BIC
[1] 418.9579
```

``` r
# Bayesian Poisson hidden Markov Model
BayesPHHMM_2states <- bayes.PHMM(y = Homicides, m = 2, chains = 3, iter = 1000, 
                                 control = list(adapt_delta = 0.99))
print(BayesPHHMM_2states, digits_summary = 3)

Inference for Stan model: PHMM.
3 chains, each with iter=1000; warmup=500; thin=1; 
post-warmup draws per chain=500, total post-warmup draws=1500.

              mean se_mean    sd     2.5%      25%      50%      75%    97.5% n_eff  Rhat
A[1,1]       0.953   0.001 0.033    0.872    0.937    0.961    0.978    0.995  1281 1.000
A[1,2]       0.047   0.001 0.033    0.005    0.022    0.039    0.063    0.128  1281 1.000
A[2,1]       0.096   0.002 0.066    0.011    0.047    0.081    0.132    0.255  1451 0.999
A[2,2]       0.904   0.002 0.066    0.745    0.868    0.919    0.953    0.989  1451 0.999
lambda[1]   29.700   0.027 0.897   27.920   29.070   29.745   30.312   31.344  1138 0.999
lambda[2]   62.823   0.053 1.922   59.038   61.594   62.758   64.078   66.619  1311 1.000
lp__      -210.670   0.062 1.517 -214.439 -211.417 -210.329 -209.543 -208.827   596 1.001

Samples were drawn using NUTS(diag_e) at Tue Feb 16 13:15:18 2021.
For each parameter, n_eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor on split chains (at 
convergence, Rhat=1).
```
