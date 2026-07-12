# Keep samples based on metadata column criteria

Keep samples based on metadata column criteria

## Usage

``` r
# S3 method for class 'nmr_dataset_family'
filter(.data, ...)
```

## Arguments

- .data:

  An
  [nmr_dataset_family](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_family.md)
  object

- ...:

  conditions, as in
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)

## Value

The same object, with the matching rows

## See also

Other subsetting functions:
[`[.nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset.md),
[`[.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset_1D.md),
[`[.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset_peak_table.md),
[`nmr_pca_outliers_filter()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_filter.md)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))

## example 1
sample_10 <- filter(dataset_1D, NMRExperiment == "10")

## example 2
# test_samples <- dataset_1D %>% filter(nmr_peak_table$metadata$external$Group == "placebo")
```
