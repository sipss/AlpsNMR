# Compute PCA residuals and score distance for each sample

Compute PCA residuals and score distance for each sample

## Usage

``` r
nmr_pca_outliers(
  nmr_dataset,
  pca_model,
  ncomp = NULL,
  quantile_critical = 0.975
)
```

## Arguments

- nmr_dataset:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object

- pca_model:

  A pca model returned by
  [nmr_pca_build_model](https://sipss.github.io/AlpsNMR/reference/nmr_pca_build_model.md)

- ncomp:

  Number of components to use. Use `NULL` for 90% of the variance

- quantile_critical:

  critical quantile

## Value

A list with:

- outlier_info: A data frame with the NMRExperiment, the Q residuals and
  T scores

- ncomp: Number of components used to compute Q and T

- Tscore_critical, QResidual_critical: Critical values, given a
  quantile, for both Q and T.

## See also

Other PCA related functions:
[`nmr_pca_build_model()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_build_model.md),
[`nmr_pca_outliers_filter()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_filter.md),
[`nmr_pca_outliers_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_plot.md),
[`nmr_pca_outliers_robust()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_robust.md),
[`nmr_pca_plots`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_plots.md)

Other outlier detection functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`nmr_pca_outliers_filter()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_filter.md),
[`nmr_pca_outliers_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_plot.md),
[`nmr_pca_outliers_robust()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_robust.md)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
model <- nmr_pca_build_model(dataset_1D)
outliers_info <- nmr_pca_outliers(dataset_1D, model)
```
