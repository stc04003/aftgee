
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/aftgee)](https://cran.r-project.org/package=aftgee)
[![packageversion](https://img.shields.io/badge/Package%20version-1.2.0-orange.svg?style=flat-square)](commits/master)
[![Travis-CI Build
Status](https://travis-ci.org/stc04003/aftgee.svg?branch=master)](https://travis-ci.org/stc04003/aftgee)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/stc04003/aftgee?branch=master&svg=true)](https://ci.appveyor.com/project/stc04003/aftgee)
[![Last-changedate](https://img.shields.io/badge/last%20change-2023--09--26-yellowgreen.svg)](/commits/master)

------------------------------------------------------------------------

## **aftgee**

------------------------------------------------------------------------

## aftgee: Accelerated Failure Time with Generalized Estimating Equations

The **aftgee** package implements recently developed inference
procedures for the accelerated failure time (AFT) models with both the
rank-based approach and the least squares approach. For the rank-based
approach, the package allows various weight choices and uses an induced
smoothing procedure that leads to much more efficient computation than
the linear programming method. With the rank-based estimator as an
initial value, the generalized estimating equation (GEE) approach is
used as an extension of the least squares approach to the multivariate
case, where the within cluster dependency is accounted for with working
covariance structures. Additional sampling weights are incorporated to
handle missing data needed as in case-cohort studies or general sampling
schemes.

## Installation

Install and load the package from CRAN using

``` r
> install.packages("aftgee")
> library(aftgee)
```

Install and load the package from GitHub using

``` r
> devtools::install_github("stc04003/aftgee")
> library(aftgee)
> packageVersion("aftgee")
```

### Online documentation

[Online document](https://www.sychiou.com/aftgee/index.html) includes:

- Package vignette on [rank-based estimators with
  `aftsrr()`](https://www.sychiou.com/aftgee/articles/aftsrr.html).
- Package vignette on [least-squares estimators with
  `aftgee()`](https://www.sychiou.com/aftgee/articles/aftgee.html).

## Reference

Chiou, S., Kang, S., and Yan, J. (2015). Rank-based estimating equations
with general weight for accelerated failure time models: An induced
smoothing approach. *Statistics in Medicine*, **34**(9): 1495–1510.

Chiou, S., Kang, S., and Yan, J. (2014). Fitting accelerated failure
time model in routine survival analysis with R package aftgee. *Journal
of Statistical Software*, **61**(11): 1–23.

Chiou, S., Kang, S., and Yan, J. (2014). Fast accelerated failure time
modeling for case-cohort data. *Statistics and Computing*, **24**(4):
559–568.

Chiou, S., Kang, S., Kim, J., and Yan, J. (2014). Marginal
semiparametric multivariate accelerated failure time model with
generalized estimating equations. *Lifetime Data Analysis*, **20**(4):
599–618.
