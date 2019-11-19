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
<<<<<<< HEAD
install.packages(c("remotes","BiocManager"))
BiocManager::install(c("impute", "MassSpecWavelet"))
=======
install.packages(c("BiocManager", "remotes"))
BiocManager::install(c("MassSpecWavelet", "impute"), update = FALSE)
>>>>>>> 1c8bae5f1de728bfc611d4f214d2cb9ea86fc373
remotes::install_github("sipss/AlpsNMR")
```

Quick start
=============

Checkout the [Introduction to AlpsNMR](https://sipss.github.io/AlpsNMR/articles/introduction-to-alpsnmr.html) vignette that shows how to import data and preprocess it using AlpsNMR.

See also the [tutorial](https://github.com/sipss/AlpsNMR/blob/master/vignettes/tutorial.pdf) with a real dataset from beginning to end, including all the steps of an untargeted metabolomics analysis. 
