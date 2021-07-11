
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/aftgee)](https://cran.r-project.org/package=aftgee)
[![packageversion](https://img.shields.io/badge/Package%20version-1.1.6-orange.svg?style=flat-square)](commits/master)
[![Travis-CI Build
Status](https://travis-ci.org/stc04003/aftgee.svg?branch=master)](https://travis-ci.org/stc04003/aftgee)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/stc04003/aftgee?branch=master&svg=true)](https://ci.appveyor.com/project/stc04003/aftgee)
[![Last-changedate](https://img.shields.io/badge/last%20change-2021--07--11-yellowgreen.svg)](/commits/master)

-----

## **aftgee**

-----

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
```

### Online documentation

[Online document](https://www.sychiou.com/aftgee/index.html).

  - Package vignette coming up.

## Examples

Here are some examples to get started:

``` r
> ## Simulate data from an AFT model with possible associated response
> datgen <- function(n = 100, tau = 0.3, dim = 2) {
+   x1 <- rbinom(dim * n, 1, 0.5)
+   x2 <- rnorm(dim * n)
+   e <- c(t(exp(MASS::mvrnorm(n = n, mu = rep(0, dim), Sigma = tau + (1 - tau) * diag(dim)))))
+   tt <- exp(2 + x1 + x2 + e)
+   cen <- runif(n, 0, 100)
+   data.frame(Time = pmin(tt, cen), status = 1 * (tt < cen),
+              x1 = x1, x2 = x2, id = rep(1:n, each = dim))
+ }
> set.seed(1); dat <- datgen(n = 100, dim = 2)
```

``` r
> library(aftgee)
> set.seed(123)
> fm <- Surv(Time, status) ~ x1 + x2
> ## Fits semiparametric AFT with rank-based approach
> fit.rk <- aftsrr(fm, id = id, data = dat, se = c("ISMB", "ZLMB"))
> summary(fit.rk)
Call:
aftsrr(formula = fm, data = dat, id = id, se = c("ISMB", "ZLMB"))

Variance Estimator: ISMB
   Estimate StdErr z.value   p.value    
x1    1.104  0.124   8.891 < 2.2e-16 ***
x2    0.963  0.075  12.870 < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Variance Estimator: ZLMB
   Estimate StdErr z.value   p.value    
x1    1.104  0.120   9.186 < 2.2e-16 ***
x2    0.963  0.085  11.305 < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
> ## Fits semiparametric AFT with GEE approach and an independent correlation structure
> fit.gee1 <- aftgee(fm, id = id, data = dat, corstr = "ind")
> summary(fit.gee1)
Call:
aftgee(formula = fm, data = dat, id = id, corstr = "ind")

AFTGEE Estimator
            Estimate StdErr z.value   p.value    
(Intercept)    3.369  0.116  28.994 < 2.2e-16 ***
x1             1.054  0.145   7.288 < 2.2e-16 ***
x2             0.958  0.081  11.890 < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
> ## Fits semiparametric AFT with GEE approach and an exchangeable correlation structure
> fit.gee2 <- aftgee(fm, id = id, data = dat, corstr = "ex")
> summary(fit.gee2)
Call:
aftgee(formula = fm, data = dat, id = id, corstr = "ex")

AFTGEE Estimator
            Estimate StdErr z.value   p.value    
(Intercept)    3.369  0.115  29.196 < 2.2e-16 ***
x1             1.055  0.146   7.244 < 2.2e-16 ***
x2             0.960  0.085  11.333 < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

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
