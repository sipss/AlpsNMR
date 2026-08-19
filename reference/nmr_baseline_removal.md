# Baseline Removal NMR

Removes the baseline on an
[nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
object, using
[baseline::baseline.als](https://rdrr.io/pkg/baseline/man/baseline.als.html).

## Usage

``` r
nmr_baseline_removal(nmr_dataset, lambda = 6, p = 0.05, maxit = 20)
```

## Arguments

- nmr_dataset:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md).

- lambda:

  2nd derivative constraint

- p:

  Weighting of positive residuals

- maxit:

  Maximum number of iterations

## Value

The same
[nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
object after baseline removal.

## Details

**Deprecated**: `nmr_baseline_removal()` will be removed in early 2027.
Use
[`nmr_baseline_estimation()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_estimation.md)
instead: it estimates the baseline without overwriting `data_1r`, and
downstream functions (e.g.
[`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md),
[`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md))
pick it up automatically when present.

## See also

[baseline::baseline.als](https://rdrr.io/pkg/baseline/man/baseline.als.html)

Other baseline removal functions:
[`nmr_baseline_estimation()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_estimation.md)

## Examples

``` r
dataset_1D <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
dataset_no_base_line <- nmr_baseline_removal(dataset_1D, lambda = 6, p = 0.01)
#> Warning: ! `nmr_baseline_removal()` is deprecated and will be removed in early 2027.
#> ℹ Use `nmr_baseline_estimation()` instead.
#> This warning is displayed once every 8 hours.
```
