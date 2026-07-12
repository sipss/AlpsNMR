# Plotting functions for PCA

Plotting functions for PCA

## Usage

``` r
nmr_pca_plot_variance(pca_model)

nmr_pca_scoreplot(nmr_dataset, pca_model, comp = seq_len(2), ...)

nmr_pca_loadingplot(pca_model, comp)
```

## Arguments

- pca_model:

  A PCA model trained with
  [nmr_pca_build_model](https://sipss.github.io/AlpsNMR/reference/nmr_pca_build_model.md)

- nmr_dataset:

  an
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object

- comp:

  Components to represent

- ...:

  Additional aesthetics passed on to
  [ggplot2::aes](https://ggplot2.tidyverse.org/reference/aes.html) (use
  bare unquoted names)

## Value

Plot of PCA

## See also

Other PCA related functions:
[`nmr_pca_build_model()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_build_model.md),
[`nmr_pca_outliers()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers.md),
[`nmr_pca_outliers_filter()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_filter.md),
[`nmr_pca_outliers_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_plot.md),
[`nmr_pca_outliers_robust()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_robust.md)

## Examples

``` r
dataset_1D <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
model <- nmr_pca_build_model(dataset_1D)
nmr_pca_plot_variance(model)

nmr_pca_scoreplot(dataset_1D, model)

nmr_pca_loadingplot(model, 1)

```
