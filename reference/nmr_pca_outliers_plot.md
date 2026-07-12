# Plot for outlier detection diagnostic

Plot for outlier detection diagnostic

## Usage

``` r
nmr_pca_outliers_plot(nmr_dataset, pca_outliers, ...)
```

## Arguments

- nmr_dataset:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object

- pca_outliers:

  The output from
  [`nmr_pca_outliers()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers.md)

- ...:

  Additional parameters passed on to
  [`ggplot2::aes()`](https://ggplot2.tidyverse.org/reference/aes.html)
  (or now deprecated to
  [`ggplot2::aes_string()`](https://ggplot2.tidyverse.org/reference/aes_.html))

## Value

A plot for the outlier detection

## See also

Other PCA related functions:
[`nmr_pca_build_model()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_build_model.md),
[`nmr_pca_outliers()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers.md),
[`nmr_pca_outliers_filter()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_filter.md),
[`nmr_pca_outliers_robust()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_robust.md),
[`nmr_pca_plots`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_plots.md)

Other outlier detection functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`nmr_pca_outliers()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers.md),
[`nmr_pca_outliers_filter()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_filter.md),
[`nmr_pca_outliers_robust()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_robust.md)

## Examples

``` r
# dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
# dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
# dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
# model <- nmr_pca_build_model(dataset_1D)
# outliers_info <- nmr_pca_outliers(dataset_1D, model)
# nmr_pca_outliers_plot(dataset_1D, outliers_info)
```
