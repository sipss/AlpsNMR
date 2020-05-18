# AlpsNMR

[![Build Status](https://github.com/sipss/AlpsNMR/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/sipss/AlpsNMR/actions/) [![codecov.io](https://codecov.io/github/sipss/AlpsNMR/coverage.svg?branch=master)](https://codecov.io/github/sipss/AlpsNMR) [![Documentation](https://img.shields.io/badge/documentation-pkgdown-informational)](https://sipss.github.io/AlpsNMR/)
[![Publication](https://img.shields.io/badge/Bioinformatics-Accepted-success)](https://doi.org/10.1093/bioinformatics/btaa022)

`AlpsNMR` is an R package that can load Bruker and JDX samples as well as
preprocess them.

It includes functions for region exclusion, normalization, peak detection & integration and
outlier detection among others. See the package vignette for details.


## Installation

AlpsNMR can be installed with the `devtools` package. For this is needed Rtools and note that it uses packages from
CRAN, from BioConductor and from git repositories:

Download Rtools for version 3.6
[[Rtools 3.6](https://cran.r-project.org/bin/windows/Rtools/Rtools35.exe)]

As you can see in Rtools website, is needed one additional step, Putting Rtools on the PATH
[[Rtools web](https://cran.r-project.org/bin/windows/Rtools/)]

The easiest way is to create Renviron file executing this command in R, (take care because if exist previous Renviron file will be erased, in this case that that file to the Renvion):

```r
writeLines('PATH="C:\\Rtool\\bin;${PATH}"', con = "~/.Renviron")
```

Now restart R, and verify that make can be found, which should show the path to your Rtools installation. (In Rstudio, ctrl+shift+F10 restart R session)

```r
Sys.which("make")
## "C:\\Rtools\\bin\\make.exe"
```
Install AlpsNMR:

```r
if (!"devtools" %in% rownames(installed.packages()))\n \tinstall.packages("devtools") 
devtools::install_github("sipss/AlpsNMR")
```

Quick start
=============

Checkout the [Introduction to AlpsNMR](https://sipss.github.io/AlpsNMR/articles/introduction-to-alpsnmr.html) vignette that shows how to import data and preprocess it using AlpsNMR. See our [publication](https://doi.org/10.1093/bioinformatics/btaa022) for further details.

See also the [tutorial](https://github.com/sipss/AlpsNMR/blob/master/vignettes/tutorial.pdf) with a real dataset from beginning to end, including all the steps of untargeted metabolomics analysis. To run the [tutorial](https://github.com/sipss/AlpsNMR/blob/master/vignettes/tutorial.pdf), you can download the MTBLS242 dataset from the public [MetaboLights repository](https://www.ebi.ac.uk/metabolights/MTBLS242), or download and unzip the contents (spectra and metadata) of this [Dropbox link](https://dl.dropboxusercontent.com/s/0snivrsd7m82yey/MTBLS242.zip?dl=0).
