# AlpsNMR

[![Build Status](https://travis-ci.org/sipss/AlpsNMR.svg?branch=master)](https://travis-ci.org/sipss/AlpsNMR) [![codecov.io](https://codecov.io/github/sipss/AlpsNMR/coverage.svg?branch=master)](https://codecov.io/github/sipss/AlpsNMR) [![Documentation](https://img.shields.io/badge/documentation-pkgdown-informational)](https://sipss.github.io/AlpsNMR/)

`AlpsNMR` is an R package that can load Bruker and JDX samples as well as
preprocess them.

It includes functions for region exclusion, normalization, peak detection & integration and
outlier detection among others. See the package vignette for details.


## Installation

AlpsNMR installation requires at least R 3.5, and some packages not yet published on CRAN or BioConductor:

    install.packages(c("BiocManager", "remotes"))
    BiocManager::install(c("impute", "MassSpecWavelet", "mixOmics"))
    remotes::install_git("https://gitlab.com/CarlBrunius/MUVR.git")
    remotes::install_github("danielcanueto/rDolphin")
    remotes::install_github("sipss/AlpsNMR")

Quick start
=============

The best way to start is by running:

    browseVignettes("AlpsNMR")

