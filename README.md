NIHSnmr
=======

`NIHSnmr` is an R package that can load Bruker and JDX samples as well as
preprocess them.

It includes functions for region exclusion, normalization, peak detection & integration and
outlier detection among others. See the package vignette for details.


Installation
=============

NIHSnmr uses the [speaq](https://cran.r-project.org/web/packages/speaq/index.html) package
that depends on the [MassSpecWavelet](http://www.bioconductor.org/packages/release/bioc/html/MassSpecWavelet.html)
and [impute](http://www.bioconductor.org/packages/release/bioc/html/impute.html) bioconductor
packages.

We can install those packages through the `BiocManager` package:

    install.packages("BiocManager")
    BiocManager::install(c("impute", "MassSpecWavelet"), version = "3.8")

In case you are on R<3.5, the installation of the speaq 2.4 package may fail due to
an update of the `Rfast` package dependency. This issue can be solved by installing
an older version of `Rfast` with the command:

    if (getRversion() < "3.5") {
      remotes::install_version("Rfast", "1.8.9")
    }

Then we can use the `remotes` package to install the rest of the required
dependencies and then build and install our package:

    remotes::install_local("NIHSnmr_2.3.0.tar.gz")

Quick start
=============

The best way to start is by running:

    browseVignettes("NIHSnmr")

