# Validate nmr_dataset_peak_table objects

Validate nmr_dataset_peak_table objects

## Usage

``` r
validate_nmr_dataset_peak_table(nmr_dataset_peak_table)
```

## Arguments

- nmr_dataset_peak_table:

  An
  [nmr_dataset_peak_table](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_peak_table.md)
  object

## Value

The
[nmr_dataset_peak_table](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_peak_table.md)
unchanged

## See also

Other class helper functions:
[`format.nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset.md),
[`format.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_1D.md),
[`format.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_peak_table.md),
[`is.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_1D.md),
[`is.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_peak_table.md),
[`new_nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/new_nmr_dataset.md),
[`new_nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/new_nmr_dataset_1D.md),
[`new_nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/new_nmr_dataset_peak_table.md),
[`print.nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset.md),
[`print.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_1D.md),
[`print.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_peak_table.md),
[`validate_nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset.md),
[`validate_nmr_dataset_family()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset_family.md)

## Examples

``` r
pt <- new_nmr_dataset_peak_table(
    peak_table = matrix(c(1, 2), nrow = 1, dimnames = list("10", c("ppm_1.4", "ppm_1.6"))),
    metadata = list(external = data.frame(NMRExperiment = "10"))
)
pt_validated <- validate_nmr_dataset_peak_table(pt)
```
