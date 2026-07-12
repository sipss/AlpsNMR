# Export data for the ASICS spectral quantification library

Exports the spectra matrix, sample names and chemical shift axis into an
ASICS Spectra object.

## Usage

``` r
to_ASICS(dataset, ...)
```

## Arguments

- dataset:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object

- ...:

  Arguments passed on to
  [`ASICS::createSpectra`](https://rdrr.io/pkg/ASICS/man/createSpectra.html)

  `norm.method`

  : Character specifying the normalisation method to use on spectra ONLY
    if the
    [`importSpectra`](https://rdrr.io/pkg/ASICS/man/importSpectra.html)
    function was not used.

  `norm.params`

  : List containing normalisation parameteres (see
    [`normaliseSpectra`](https://rdrr.io/pkg/ASICS/man/normaliseSpectra.html)
    for details) ONLY if the
    [`importSpectra`](https://rdrr.io/pkg/ASICS/man/importSpectra.html)
    function was not used.

## Value

An [ASICS::Spectra](https://rdrr.io/pkg/ASICS/man/Spectra-class.html)
object

## Examples

``` r
if (requireNamespace("ASICS", quietly=TRUE)) {
  nsamp <- 3
  npoints <- 300
  metadata <- list(external = data.frame(
    NMRExperiment = paste0("Sample", seq_len(nsamp))
  ))
  dataset <- new_nmr_dataset_1D(
    ppm_axis = seq(from = 0.2, to = 10, length.out = npoints),
    data_1r = matrix(runif(nsamp * npoints), nrow = nsamp, ncol = npoints),
    metadata = metadata
  )
  forAsics <- to_ASICS(dataset)
  #ASICS::ASICS(forAsics)
}
```
