# Build a PCA on for an nmr_dataset

This function builds a PCA model with all the NMR spectra. Regions with
zero values (excluded regions) or near-zero variance regions are
automatically excluded from the analysis.

## Usage

``` r
nmr_pca_build_model(
  nmr_dataset,
  ncomp = NULL,
  center = TRUE,
  scale = FALSE,
  ...
)

# S3 method for class 'nmr_dataset_1D'
nmr_pca_build_model(
  nmr_dataset,
  ncomp = NULL,
  center = TRUE,
  scale = FALSE,
  ...
)
```

## Arguments

- nmr_dataset:

  a
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object

- ncomp:

  Integer, if data is complete `ncomp` decides the number of components
  and associated eigenvalues to display from the `pcasvd` algorithm and
  if the data has missing values, `ncomp` gives the number of components
  to keep to perform the reconstitution of the data using the NIPALS
  algorithm. If `NULL`, function sets `ncomp = min(nrow(X), ncol(X))`

- center:

  (Default=TRUE) Logical, whether the variables should be shifted to be
  zero centered. Only set to FALSE if data have already been centered.
  Alternatively, a vector of length equal the number of columns of `X`
  can be supplied. The value is passed to
  [`scale`](https://rdrr.io/r/base/scale.html). If the data contain
  missing values, columns should be centered for reliable results.

- scale:

  (Default=FALSE) Logical indicating whether the variables should be
  scaled to have unit variance before the analysis takes place. The
  default is `FALSE` for consistency with `prcomp` function, but in
  general scaling is advisable. Alternatively, a vector of length equal
  the number of columns of `X` can be supplied. The value is passed to
  [`scale`](https://rdrr.io/r/base/scale.html).

- ...:

  Additional arguments passed on to
  [mixOmics::pca](https://rdrr.io/pkg/mixOmics/man/pca.html)

## Value

A PCA model as given by
[mixOmics::pca](https://rdrr.io/pkg/mixOmics/man/pca.html) with two
additional attributes:

- `nmr_data_axis` containing the full ppm axis

- `nmr_included` with the data points included in the model These
  attributes are used internally by AlpsNMR to create loading plots

## See also

Other PCA related functions:
[`nmr_pca_outliers()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers.md),
[`nmr_pca_outliers_filter()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_filter.md),
[`nmr_pca_outliers_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_plot.md),
[`nmr_pca_outliers_robust()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_robust.md),
[`nmr_pca_plots`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_plots.md)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
model <- nmr_pca_build_model(dataset_1D)

dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
model <- nmr_pca_build_model(dataset_1D)
```
