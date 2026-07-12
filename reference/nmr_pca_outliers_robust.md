# Outlier detection through robust PCA

Outlier detection through robust PCA

## Usage

``` r
nmr_pca_outliers_robust(nmr_dataset, ncomp = 5)
```

## Arguments

- nmr_dataset:

  An nmr_dataset_1D object

- ncomp:

  Number of rPCA components to use

  We have observed that the statistical test used as a threshold for
  outlier detection usually flags as outliers too many samples, due
  possibly to a lack of gaussianity

  As a workaround, a heuristic method has been implemented: We know that
  in the Q residuals vs T scores plot from
  [`nmr_pca_outliers_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_plot.md)
  outliers are on the right or on the top of the plot, and quite
  separated from non-outlier samples.

  To determine the critical value, both for Q and T, we find the biggest
  gap between samples in the plot and use as critical value the center
  of the gap.

  This approach seems to work well when there are outliers, but it fails
  when there isn't any outlier. For that case, the gap would be placed
  anywhere and that is not desirable as many samples would be
  incorrectly flagged. The second assumption that we use is that no more
  than 10\\ the samples may pass our critical value. If more than 10\\
  pass the critical value, then we assume that our heuristics are not
  reasonable and we don't set any critical limit.

## Value

A list similar to
[nmr_pca_outliers](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers.md)

## See also

Other PCA related functions:
[`nmr_pca_build_model()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_build_model.md),
[`nmr_pca_outliers()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers.md),
[`nmr_pca_outliers_filter()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_filter.md),
[`nmr_pca_outliers_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_plot.md),
[`nmr_pca_plots`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_plots.md)

Other outlier detection functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`nmr_pca_outliers()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers.md),
[`nmr_pca_outliers_filter()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_filter.md),
[`nmr_pca_outliers_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_plot.md)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
outliers_info <- nmr_pca_outliers_robust(dataset_1D)
```
