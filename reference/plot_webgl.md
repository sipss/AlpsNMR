# Plot a dataset into a HTML file

Uses WebGL for performance

## Usage

``` r
plot_webgl(nmr_dataset, html_filename, overwrite = NULL, ...)
```

## Arguments

- nmr_dataset:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)

- html_filename:

  The output HTML filename to be created

- overwrite:

  Overwrite the lib/ directory (use `NULL` to prompt the user)

- ...:

  Arguments passed on to
  [`plot.nmr_dataset_1D`](https://sipss.github.io/AlpsNMR/reference/plot.nmr_dataset_1D.md)

  `x`

  : a
    [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
    object

  `chemshift_range`

  : range of the chemical shifts to be included. Can be of length 3 to
    include the resolution in the third element (e.g.
    `c(0.2, 0.8, 0.005)`)

  `NMRExperiment`

  : A character vector with the NMRExperiments to include. Use "all" to
    include all experiments.

  `quantile_plot`

  : If `TRUE` plot the 10\\ If two numbers between 0 and 1 are given
    then a custom percentile can be plotted

  `quantile_colors`

  : A vector with the colors for each of the quantiles

  `interactive`

  : if `TRUE` return an interactive plotly plot, otherwise return a
    ggplot one.

## Value

the html filename created

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
# dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
# dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
# html_plot <- plot_webgl(dataset_1D, "html_plot.html")
```
