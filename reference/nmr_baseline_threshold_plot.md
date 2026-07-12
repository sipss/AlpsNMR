# Plot the baseline thresholds

If you have a lot of samples you can make the plot window bigger (or use
"````  ```{r fig.height=10, fig.width=10} ````" in notebooks), or choose
some NMRExperiments.

## Usage

``` r
nmr_baseline_threshold_plot(
  nmr_dataset,
  thresholds,
  NMRExperiment = "all",
  chemshift_range = c(9.5, 10),
  ...
)
```

## Arguments

- nmr_dataset:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object

- thresholds:

  A named vector. The values are baseline thresholds. The names are
  NMRExperiments.

- NMRExperiment:

  The NMRExperiments to plot (Use `"all"` to plot all of them)

- chemshift_range:

  The range to plot, as a first check use the `range_without_peaks` from
  [nmr_baseline_threshold](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md)

- ...:

  arguments passed to
  [ggplot2::aes](https://ggplot2.tidyverse.org/reference/aes.html) (or
  to
  [ggplot2::aes_string](https://ggplot2.tidyverse.org/reference/aes_.html),
  being deprecated).

## Value

A plot.

## Examples

``` r

ppm_axis <- seq(from = 0, to = 10, length.out = 1000)
data_1r <- matrix(runif(1000, 0, 10), nrow = 1) + 100
dataset_1D <- new_nmr_dataset_1D(
    ppm_axis = ppm_axis,
    data_1r = data_1r,
    metadata = list(external=data.frame(NMRExperiment = "10"))
)
bl_threshold <- nmr_baseline_threshold(dataset_1D, range_without_peaks = c(9.5,10))
baselineThresh <- nmr_baseline_threshold(dataset_1D)
nmr_baseline_threshold_plot(dataset_1D, bl_threshold)
```
