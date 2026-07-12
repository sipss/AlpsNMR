# Rephase 1D NMR data

Use phasing algorithms to rephase data in the spectral domain.

This function may improve autophasing processing from instrument
vendors. It wraps the
[`NMRphasing::NMRphasing()`](https://rdrr.io/pkg/NMRphasing/man/NMRphasing.html)
function, to automatically rephase spectra, allowing you to choose from
a number of algorithms of which `NLS`, `MPC_DANM` and `SPC_DANM` are the
most recent.

Rephasing should happen before any spectra interpolation.

Please use the `all_components = TRUE` when calling
[`nmr_read_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_samples.md)
in order to load the complex spectra and fix NMR phasing correctly.

## Usage

``` r
nmr_autophase(
  dataset,
  method = c("NLS", "MPC_DANM", "MPC_EMP", "SPC_DANM", "SPC_EMP", "SPC_AAM", "SPC_DSM"),
  withBC = FALSE,
  ...
)
```

## Arguments

- dataset:

  An
  [nmr_dataset](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md)
  object

- method:

  The autophasing method. See
  [`NMRphasing::NMRphasing()`](https://rdrr.io/pkg/NMRphasing/man/NMRphasing.html)
  for details.

- withBC:

  [`NMRphasing::NMRphasing`](https://rdrr.io/pkg/NMRphasing/man/NMRphasing.html)
  may perform a baseline correction using modified polynomial fitting.
  By default AlpsNMR offers other baseline estimation methods and better
  visualization of its effect, so AlpsNMR by default disables the
  baseline correction offered by NMRphasing.

- ...:

  Other parameters passed on to
  [`NMRphasing::NMRphasing()`](https://rdrr.io/pkg/NMRphasing/man/NMRphasing.html).

## Value

A (hopefully better phased)
[nmr_dataset](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md)
object, with updated real and imaginary parts.

## Examples

``` r
if (requireNamespace("NMRphasing", quietly=TRUE)) {
  # Helpers to create a dataset:
  lorentzian <- function(x, x0, gamma, A) {
    A * (1 / (pi * gamma)) * ((gamma^2) / ((x - x0)^2 + gamma^2))
  }
  x <- seq(from=1, to=2, length.out = 300)
  y <- lorentzian(x, 1.3, 0.01, 1) + lorentzian(x, 1.6, 0.01, 1)
  dataset <- new_nmr_dataset(
    metadata = list(external = data.frame(NMRExperiment = "10")),
    data_fields = list(data_1r = list(y)),
    axis = list(list(x))
  )
  # Autophase, interpolate and plot:
  dataset <- nmr_autophase(dataset, method = "NLS")
  dataset <- nmr_interpolate_1D(dataset, axis = c(min = 1, max = 2, by = 0.01))
  plot(dataset)
}
#> Warning: ! nmr_autophase() performs better with access to the whole complex NMR spectra
#> ℹ Please read the dataset using `all_components=TRUE`.
#> ℹ See `AlpsNMR::nmr_autophase()` for a full example
```
