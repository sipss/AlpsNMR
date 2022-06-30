# AlpsNMR <img src='man/figures/AlpsNMRlogo.png' align="right" width="120" height="139" />

[![Build Status](https://github.com/sipss/AlpsNMR/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/sipss/AlpsNMR/actions/)
[![codecov.io](https://codecov.io/github/sipss/AlpsNMR/coverage.svg?branch=master)](https://codecov.io/github/sipss/AlpsNMR)
[![Bioc Status](https://bioconductor.org/shields/build/devel/bioc/AlpsNMR.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/AlpsNMR/)
[![Documentation](https://img.shields.io/badge/documentation-pkgdown-informational)](https://sipss.github.io/AlpsNMR/)
[![Publication](https://img.shields.io/badge/Bioinformatics-Accepted-success)](https://doi.org/10.1093/bioinformatics/btaa022)

`AlpsNMR` is an R package that can load Bruker and JDX samples as well as
preprocess them.

It includes functions for region exclusion, normalization, peak detection & integration and
outlier detection among others. See the package vignette for details.

## Installation

### Latest release

```r
if (!"BiocManager" %in% rownames(installed.packages()))  
    install.packages("BiocManager")
BiocManager::install("AlpsNMR")
```

### Development version

```r
if (!"remotes" %in% rownames(installed.packages()))  
    install.packages("remotes")
remotes::install_github("sipss/AlpsNMR")
```


Quick start
=============

Checkout the [Introduction to AlpsNMR](https://sipss.github.io/AlpsNMR/articles/Vig01-introduction-to-alpsnmr.html) vignette that shows how to import data and preprocess it using `AlpsNMR`. See our [publication](https://doi.org/10.1093/bioinformatics/btaa022) for further details.

See also the [tutorial](https://github.com/sipss/AlpsNMRWorkflow) with a real dataset from beginning to end, including all the steps of untargeted metabolomics analysis.
