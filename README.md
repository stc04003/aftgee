**aftgee**
----------

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/aftgee)](https://cran.r-project.org/package=aftgee) [![packageversion](https://img.shields.io/badge/Package%20version-1.1.3-orange.svg?style=flat-square)](commits/master) [![Last-changedate](https://img.shields.io/badge/last%20change-2018--06--30-yellowgreen.svg)](/commits/master)

------------------------------------------------------------------------

aftgee: Accerated Failure Time with Generalized Estimating Equation
-------------------------------------------------------------------

The **aftgee** package implements recently developed inference procedures for the accelerated failure time (AFT) models with both the rank-based approach and the least squares approach. For the rank-based approach, the package allows various weight choices and uses an induced smoothing procedure that leads to much more efficient computation than the linear programming method. With the rank-based estimator as an initial value, the generalized estimating equation (GEE) approach is used as an extension of the least squares approach to the multivariate case, where the within cluster dependency is accountered for with working covariance structures. Additional sampling weights are incorporated to handle missing data needed as in case-cohort studies or general sampling schemes.

Installation
------------

Install and load the package from CRAN using

``` r
install.packages("aftgee")
library(aftgee)
```

Install and load the package from GitHub using

``` r
devtools::install_github("stc04003/aftgee")
library(aftgee)
```

Examples
--------

Here are some examples to get started:

``` r
library(survival)
data(kidney)

library(aftgee)
fit.rk <- aftsrr(Surv(time, status) ~ age + sex, id = id, data = kidney, se = c("ISMB", "ZLMB"))
summary(fit.rk)
#> Call:
#> aftsrr(formula = Surv(time, status) ~ age + sex, data = kidney, 
#>     id = id, se = c("ISMB", "ZLMB"))
#> 
#> Variance Estimator: ISMB
#>     Estimate StdErr z.value p.value   
#> age   -0.001  0.016  -0.076   0.939   
#> sex    1.522  0.521   2.920   0.003 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Variance Estimator: ZLMB
#>     Estimate StdErr z.value p.value   
#> age   -0.001  0.017  -0.074   0.941   
#> sex    1.522  0.538   2.828   0.005 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
fit.ge <- aftgee(Surv(time, status) ~ age + sex, id = id, data = kidney)
summary(fit.ge)
#> Call:
#> aftgee(formula = Surv(time, status) ~ age + sex, data = kidney, 
#>     id = id)
#> 
#> AFTGEE Estimator
#>             Estimate StdErr z.value p.value    
#> (Intercept)    2.071  0.757   2.734   0.006 ** 
#> age           -0.005  0.007  -0.702   0.482    
#> sex            1.374  0.366   3.751  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
