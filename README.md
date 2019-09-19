# AlpsNMR


`AlpsNMR` is an R package that can load Bruker and JDX samples as well as
preprocess them.

It includes functions for region exclusion, normalization, peak detection & integration and
outlier detection among others. See the package vignette for details.


## Installation


If the installation failed, the reason is probably a missing dependency. If so,
please check the console for warning messages and install these dependecies
manually.
AlpsNMR uses the [speaq](https://cran.r-project.org/web/packages/speaq/index.html) package
that depends on the [MassSpecWavelet](http://www.bioconductor.org/packages/release/bioc/html/MassSpecWavelet.html)
and [impute](http://www.bioconductor.org/packages/release/bioc/html/impute.html) bioconductor
packages.

### Bioconductor dependencies

#### On R-release (currently R-3.5)

We can install those packages through the `BiocManager` package:

    install.packages("BiocManager")
    BiocManager::install(c("impute", "MassSpecWavelet"), version = "3.8")
    BiocManager::install("mixOmics")



#### On older R versions

    source("https://bioconductor.org/biocLite.R")
    BiocInstaller::biocLite(c("impute", "MassSpecWavelet"))

### CRAN dependencies

#### If R < 3.5

In case you are on R<3.5, the installation of the speaq 2.4 package may fail due to
an update of the `Rfast` package dependency. This issue can be solved by installing
an older version of `Rfast` with the command:

    install.packages("remotes")
    if (getRversion() < "3.5") {
      remotes::install_version("Rfast", "1.8.9")
    }

### GitLab dependencies

MUVR package can be installed with this comand:

    library(devtools)
    install_git("https://gitlab.com/CarlBrunius/MUVR.git")


### Github dependencies

rDolphin package can be installed with this comand:

    devtools::install_github("danielcanueto/rDolphin")


### AlpsNMR installation

Then we can use the `remotes` package to install the rest of the required
dependencies and then build and install our package:

    remotes::install_local("AlpsNMR_2.3.4.tar.gz")

Quick start
=============

The best way to start is by running:

    browseVignettes("AlpsNMR")

