# AlpsNMR

[![Build Status](https://travis-ci.org/sipss/AlpsNMR.svg?branch=master)](https://travis-ci.org/sipss/AlpsNMR) [![codecov.io](https://codecov.io/github/sipss/AlpsNMR/coverage.svg?branch=master)](https://codecov.io/github/sipss/AlpsNMR) [![Documentation](https://img.shields.io/badge/documentation-pkgdown-informational)](https://sipss.github.io/AlpsNMR/)

`AlpsNMR` is an R package that can load Bruker and JDX samples as well as
preprocess them.

It includes functions for region exclusion, normalization, peak detection & integration and
outlier detection among others. See the package vignette for details.


## Installation

AlpsNMR can be installed with the `remotes` package. Note that it uses packages from
CRAN, from BioConductor and from git repositories:

```r
install.packages("remotes")
remotes::install_github("sipss/AlpsNMR")
```

Quick start
=============

Checkout the [Introduction to AlpsNMR](https://sipss.github.io/AlpsNMR/articles/introduction-to-alpsnmr.html) vignette that shows how to import data and preprocess it using AlpsNMR.

Also the [Turorial](https://sipss.github.io/AlpsNMR/articles/tutorial.pdf)