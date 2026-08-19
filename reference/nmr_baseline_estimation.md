# Estimate the baseline on an [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md) object, using [`psalsa()`](https://sipss.github.io/AlpsNMR/reference/psalsa.md)

Estimates the baseline of every sample in `nmr_dataset` with the PSALSA
algorithm (see
[`psalsa()`](https://sipss.github.io/AlpsNMR/reference/psalsa.md)) and
stores it in the `data_1r_baseline` element, leaving `data_1r` itself
untouched. Several other functions
([`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md),
[`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md),
[`nmr_normalize()`](https://sipss.github.io/AlpsNMR/reference/nmr_normalize.md))
pick up `data_1r_baseline` automatically when it is present.

## Usage

``` r
nmr_baseline_estimation(
  nmr_dataset,
  lambda = "auto",
  p = "auto",
  k = "auto",
  maxit = "auto"
)
```

## Arguments

- nmr_dataset:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md).

- lambda:

  Smoothing parameter, or `"auto"` to pick it with
  [`tune_psalsa()`](https://sipss.github.io/AlpsNMR/reference/tune_psalsa.md).
  See [`psalsa()`](https://sipss.github.io/AlpsNMR/reference/psalsa.md).

- p:

  Asymmetry parameter, or `"auto"` to pick it with
  [`tune_psalsa()`](https://sipss.github.io/AlpsNMR/reference/tune_psalsa.md).
  See [`psalsa()`](https://sipss.github.io/AlpsNMR/reference/psalsa.md).

- k:

  Peak height parameter, or `"auto"` to pick it with
  [`tune_psalsa()`](https://sipss.github.io/AlpsNMR/reference/tune_psalsa.md).
  See [`psalsa()`](https://sipss.github.io/AlpsNMR/reference/psalsa.md).

- maxit:

  Maximum number of iterations, or `"auto"` to use
  [`psalsa()`](https://sipss.github.io/AlpsNMR/reference/psalsa.md)'s
  own default.

## Value

The same
[nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
object with the `data_1r_baseline` element.

## Details

`lambda`, `p` and `k` each default to `"auto"`. Whenever any of them is
`"auto"`,
[`tune_psalsa()`](https://sipss.github.io/AlpsNMR/reference/tune_psalsa.md)
is run once on every sample in `nmr_dataset` (pooled together, as it
would be for a single call to
[`tune_psalsa()`](https://sipss.github.io/AlpsNMR/reference/tune_psalsa.md)
with a list of spectra) to pick values for every `"auto"` parameter; a
parameter given as an explicit number instead bypasses tuning for that
parameter and is passed to
[`psalsa()`](https://sipss.github.io/AlpsNMR/reference/psalsa.md) as-is.
The same `lambda`/`p`/`k` are then used for every sample. `maxit` also
defaults to `"auto"`, meaning
[`psalsa()`](https://sipss.github.io/AlpsNMR/reference/psalsa.md)'s own
default is used, since `maxit` is not tuned by
[`tune_psalsa()`](https://sipss.github.io/AlpsNMR/reference/tune_psalsa.md).

## See also

[`psalsa()`](https://sipss.github.io/AlpsNMR/reference/psalsa.md),
[`tune_psalsa()`](https://sipss.github.io/AlpsNMR/reference/tune_psalsa.md)

Other baseline removal functions:
[`nmr_baseline_removal()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_removal.md)

## Examples

``` r
dataset_1D <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
dataset_1D <- nmr_baseline_estimation(dataset_1D)
```
