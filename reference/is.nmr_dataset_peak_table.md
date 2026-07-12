# Object is of [nmr_dataset_peak_table](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_peak_table.md) class

Object is of
[nmr_dataset_peak_table](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_peak_table.md)
class

## Usage

``` r
is.nmr_dataset_peak_table(x)
```

## Arguments

- x:

  an
  [nmr_dataset_peak_table](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_peak_table.md)
  object

## Value

`TRUE` if the object is an `nmr_dataset_peak_table`, `FALSE` otherwise

## See also

Other class helper functions:
[`format.nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset.md),
[`format.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_1D.md),
[`format.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_peak_table.md),
[`is.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_1D.md),
[`new_nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/new_nmr_dataset.md),
[`new_nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/new_nmr_dataset_1D.md),
[`new_nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/new_nmr_dataset_peak_table.md),
[`print.nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset.md),
[`print.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_1D.md),
[`print.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_peak_table.md),
[`validate_nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset.md),
[`validate_nmr_dataset_family()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset_family.md),
[`validate_nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset_peak_table.md)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
meta <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
metadata <- readxl::read_excel(meta, sheet = 1)
dataset_1D <- nmr_meta_add(dataset_1D, metadata = metadata, by = "NMRExperiment")
metadata <- list(external = dataset_1D[["metadata"]][["external"]])
peak_table <- nmr_data(dataset_1D)
new <- new_nmr_dataset_peak_table(peak_table, metadata)
is(new)
#> [1] "nmr_dataset_peak_table"
```
