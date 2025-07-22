
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ppls: Penalized Partial Least Squares for Structured and Functional Data

<!-- badges: start -->

<a href="https://cran.r-project.org/package=ppls"
class="ppls-release"><img
src="https://www.r-pkg.org/badges/version/ppls" alt="CRAN Status" /></a>
<!-- badges: end -->

## Introduction

This R package implements a flexible and powerful framework for
**penalized Partial Least Squares (PPLS)**, including:

- B-spline basis transformations of the input data,
- Construction of penalty matrices,
- Cross-validation for hyperparameter tuning (`lambda`, `ncomp`),
- Visualization of component effects,
- Jackknife-based inference (covariance, t-tests),
- Support for block-wise variable selection and kernel PPLS.

It is particularly suited for **high-dimensional, structured, and
nonlinear regression problems**, such as functional data or spectral
data.

------------------------------------------------------------------------

## Installation

We recommend to install the CRAN version for a stable version of `ppls`.

``` r
install.packages("ppls")
```

Alternatively, install the development version from GitHub with:

``` r
# install.packages("pak") ## If necessary, install "pak" beforehand
pak::pak("vguillemot/ppls")
```

------------------------------------------------------------------------

## Example: Fit a Penalized PLS Model with Splines

The following is an example of how use `ppls` on the `cookie` dataset
included in the package:

``` r
library(ppls)

# Load cookie data
data(cookie)
X <- as.matrix(cookie$NIR)
y <- cookie$constituents$fat

# Fit a 2-component kernel PPLS 
fit.kpls <- penalized.pls(
  X = X, y = y,
  kernel = TRUE, 
  ncomp = 2
)

# Predict on train data
yhat <- new.penalized.pls(ppls = fit.kpls, Xtest = X, ytest = y)$ypred

# Plot predicted vs observed
plot(yhat[,1], y, xlab = "Fitted", ylab = "Observed", pch = 16, asp = 1)
abline(0, 1, col = "blue", lty = 2)
```

------------------------------------------------------------------------

## References

> Kraemer, N., Boulesteix, A.-L., & Tutz, G. (2008).  
> *Penalized Partial Least Squares with Applications to B-Spline
> Transformations and Functional Data*.  
> Chemometrics and Intelligent Laboratory Systems, 94(1), 60–69.  
> <https://doi.org/10.1016/j.chemolab.2008.06.009>
