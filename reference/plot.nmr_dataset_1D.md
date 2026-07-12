# Plot an nmr_dataset_1D

Plot an nmr_dataset_1D

## Usage

``` r
# S3 method for class 'nmr_dataset_1D'
plot(
  x,
  NMRExperiment = NULL,
  chemshift_range = NULL,
  interactive = FALSE,
  quantile_plot = NULL,
  quantile_colors = NULL,
  ...
)
```

## Arguments

- x:

  a
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object

- NMRExperiment:

  A character vector with the NMRExperiments to include. Use "all" to
  include all experiments.

- chemshift_range:

  range of the chemical shifts to be included. Can be of length 3 to
  include the resolution in the third element (e.g.
  `c(0.2, 0.8, 0.005)`)

- interactive:

  if `TRUE` return an interactive plotly plot, otherwise return a ggplot
  one.

- quantile_plot:

  If `TRUE` plot the 10\\ If two numbers between 0 and 1 are given then
  a custom percentile can be plotted

- quantile_colors:

  A vector with the colors for each of the quantiles

- ...:

  arguments passed to
  [ggplot2::aes](https://ggplot2.tidyverse.org/reference/aes.html) (or
  to
  [ggplot2::aes_string](https://ggplot2.tidyverse.org/reference/aes_.html),
  being deprecated).

## Value

The plot

## See also

Other plotting functions:
[`plot_interactive()`](https://sipss.github.io/AlpsNMR/reference/plot_interactive.md)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
# dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
# dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
# plot(dataset_1D)
```
