# AlpsNMR

[![Build Status](https://github.com/sipss/AlpsNMR/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/sipss/AlpsNMR/workflows/R-CMD-check/) [![codecov.io](https://codecov.io/github/sipss/AlpsNMR/coverage.svg?branch=master)](https://codecov.io/github/sipss/AlpsNMR) [![Documentation](https://img.shields.io/badge/documentation-pkgdown-informational)](https://sipss.github.io/AlpsNMR/)

`AlpsNMR` is an R package that can load Bruker and JDX samples as well as
preprocess them.

It includes functions for region exclusion, normalization, peak detection & integration and
outlier detection among others. See the package vignette for details.


## Installation

AlpsNMR can be installed with the `remotes` package. Note that it uses packages from
CRAN, from BioConductor and from git repositories:

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install(c("MassSpecWavelet", "impute"), update = FALSE)
remotes::install_github("sipss/AlpsNMR")
```

Quick start
=============

Checkout the [Introduction to AlpsNMR](https://sipss.github.io/AlpsNMR/articles/introduction-to-alpsnmr.html) vignette that shows how to import data and preprocess it using AlpsNMR.

See also the [tutorial](https://github.com/sipss/AlpsNMR/blob/master/vignettes/tutorial.pdf) with a real dataset from beginning to end, including all the steps of untargeted metabolomics analysis. To run the [tutorial](https://github.com/sipss/AlpsNMR/blob/master/vignettes/tutorial.pdf), you can download the MTBLS242 dataset from the public [MetaboLights repository](https://www.ebi.ac.uk/metabolights/MTBLS242), or download and unzip the contents (spectra and metadata) of this [Dropbox link](https://dl.dropboxusercontent.com/s/0snivrsd7m82yey/MTBLS242.zip?dl=0).

