# Extract parts of an nmr_dataset_1D

Extract parts of an nmr_dataset_1D

## Usage

``` r
# S3 method for class 'nmr_dataset_1D'
x[i]
```

## Arguments

- x:

  an
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object

- i:

  indices of the samples to keep

## Value

an nmr_dataset_1D with the extracted samples

## See also

Other subsetting functions:
[`[.nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset.md),
[`[.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset_peak_table.md),
[`filter.nmr_dataset_family()`](https://sipss.github.io/AlpsNMR/reference/filter.nmr_dataset_family.md),
[`nmr_pca_outliers_filter()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_filter.md)

Other nmr_dataset_1D functions:
[`format.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_1D.md),
[`get_integration_with_metadata()`](https://sipss.github.io/AlpsNMR/reference/get_integration_with_metadata.md),
[`is.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_1D.md),
[`nmr_integrate_peak_positions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_peak_positions.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md),
[`nmr_meta_add()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_add.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_meta_get()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get.md),
[`nmr_meta_get_column()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get_column.md),
[`nmr_ppm_resolution()`](https://sipss.github.io/AlpsNMR/reference/nmr_ppm_resolution.md),
[`print.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_1D.md)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
dataset_1D[0]
#> An nmr_dataset_1D (0 samples) 
```
