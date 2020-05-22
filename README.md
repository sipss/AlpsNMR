# AlpsNMR

[![Build Status](https://github.com/sipss/AlpsNMR/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/sipss/AlpsNMR/actions/) [![codecov.io](https://codecov.io/github/sipss/AlpsNMR/coverage.svg?branch=master)](https://codecov.io/github/sipss/AlpsNMR) [![Documentation](https://img.shields.io/badge/documentation-pkgdown-informational)](https://sipss.github.io/AlpsNMR/)
[![Publication](https://img.shields.io/badge/Bioinformatics-Accepted-success)](https://doi.org/10.1093/bioinformatics/btaa022)

`AlpsNMR` is an R package that can load Bruker and JDX samples. It provides automated and efficient signal processing for untargeted 
 NMR metabolomics.

It is able to interpolate the samples, detect outliers, exclude regions, normalize, detect peaks, align the spectra, integrate peaks, manage metadata and visualize the spectra. See the package vignette for details.


## Installation

`AlpsNMR` can be installed with the `devtools` package and note that it utilizes packages from CRAN, from BioConductor, and git repositories. Therefore, `Rtools` is also required.

Probably you already have `Rtools` installed. Then, run this code to install `AlspNMR`:

```r
if (!"BiocManager" %in% rownames(installed.packages()))  
    install.packages("BiocManager")
BiocManager::install(c("MassSpecWavelet", "impute"), update = FALSE)  
if (!"devtools" %in% rownames(installed.packages()))  
    install.packages("devtools")
devtools::install_github("sipss/AlpsNMR")
```


If you need to install `Rtools`:

1. Download `Rtools` for version 3.6
[[Rtools 3.6](https://cran.r-project.org/bin/windows/Rtools/Rtools35.exe)]

2. The additional step of adding `Rtools` to the [[PATH](https://cran.r-project.org/bin/windows/Rtools/)] is recommended to ensure the proper installation of the package.

The easiest way to do this is to create an `Renviron` file executing this command in R (note that if this file already exists, it will be replaced):

```r
writeLines('PATH="C:\\Rtool\\bin;${PATH}"', con = "~/.Renviron")
```

3. Now it is time to restart R, and verify that `make` can be found. Ths command should show the path to your `Rtools` installation.

```r
Sys.which("make")
## "C:\\Rtools\\bin\\make.exe"
```

4. Install AlpsNMR:

```r
if (!"BiocManager" %in% rownames(installed.packages()))  
    install.packages("BiocManager")  
BiocManager::install(c("MassSpecWavelet", "impute"), update = FALSE)  
if (!"devtools" %in% rownames(installed.packages()))  
    install.packages("devtools")  
devtools::install_github("sipss/AlpsNMR")
```

Quick start
=============

Checkout the [Introduction to AlpsNMR](https://sipss.github.io/AlpsNMR/articles/introduction-to-alpsnmr.html) vignette that shows how to import data and preprocess it using AlpsNMR. See our [publication](https://doi.org/10.1093/bioinformatics/btaa022) for further details.

See also the [tutorial](https://github.com/sipss/AlpsNMR/blob/master/vignettes/tutorial.pdf) with a real dataset from beginning to end, including all the steps of untargeted metabolomics analysis. To run the [tutorial](https://github.com/sipss/AlpsNMR/blob/master/vignettes/tutorial.pdf), you can download the MTBLS242 dataset from the public [MetaboLights repository](https://www.ebi.ac.uk/metabolights/MTBLS242), or download and unzip the contents (spectra and metadata) of this [Dropbox link](https://dl.dropboxusercontent.com/s/0snivrsd7m82yey/MTBLS242.zip?dl=0).
