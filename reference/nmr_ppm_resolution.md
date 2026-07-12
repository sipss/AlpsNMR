# PPM resolution of the spectra

The function gets the ppm resolution of the dataset using the median of
the difference of data points.

## Usage

``` r
nmr_ppm_resolution(nmr_dataset)

# S3 method for class 'nmr_dataset'
nmr_ppm_resolution(nmr_dataset)

# S3 method for class 'nmr_dataset_1D'
nmr_ppm_resolution(nmr_dataset)
```

## Arguments

- nmr_dataset:

  An object containing NMR samples

## Value

Numeric (the ppm resolution, measured in ppms)

## See also

Other nmr_dataset_1D functions:
[`[.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset_1D.md),
[`format.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_1D.md),
[`get_integration_with_metadata()`](https://sipss.github.io/AlpsNMR/reference/get_integration_with_metadata.md),
[`is.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_1D.md),
[`nmr_integrate_peak_positions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_peak_positions.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md),
[`nmr_meta_add()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_add.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_meta_get()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get.md),
[`nmr_meta_get_column()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get_column.md),
[`print.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_1D.md)

## Examples

``` r
nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
nmr_ppm_resolution(nmr_dataset)
#> [1] 0.02
message("the ppm resolution of this dataset is ", nmr_ppm_resolution(nmr_dataset), " ppm")
#> the ppm resolution of this dataset is 0.02 ppm

nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
nmr_ppm_resolution(nmr_dataset)
#> [1] 0.02
message("the ppm resolution of this dataset is ", nmr_ppm_resolution(nmr_dataset), " ppm")
#> the ppm resolution of this dataset is 0.02 ppm

nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
nmr_ppm_resolution(nmr_dataset)
#> [1] 0.02
message("the ppm resolution of this dataset is ", nmr_ppm_resolution(nmr_dataset), " ppm")
#> the ppm resolution of this dataset is 0.02 ppm
```
