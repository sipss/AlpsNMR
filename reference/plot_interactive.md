# Plots in WebGL

Plots in WebGL

## Usage

``` r
plot_interactive(plt, html_filename, overwrite = NULL)
```

## Arguments

- plt:

  A plot created with plotly or ggplot2

- html_filename:

  The file name where the plot will be saved

- overwrite:

  Overwrite the lib/ directory (use `NULL` to prompt the user)

## Value

The html_filename

## See also

Other plotting functions:
[`plot.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/plot.nmr_dataset_1D.md)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
# plot <- plot(dataset_1D)
# html_plot_interactive <- plot_interactive(plot, "html_plot_interactive.html")
```
