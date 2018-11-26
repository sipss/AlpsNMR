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

Then we can use the `remotes` package to install the rest of the required
dependencies and then build and install our package:

    remotes::install_local("NIHSnmr_2.1.0.tar.gz")

Quick start
=============

The best way to start is by running:

    browseVignettes("NIHSnmr")

