# Plot clustering results

Plot clustering results

## Usage

``` r
nmr_peak_clustering_plot(
  dataset,
  peak_list_clustered,
  NMRExperiments,
  chemshift_range,
  baselineThresh = NULL
)
```

## Arguments

- dataset:

  The
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object

- peak_list_clustered:

  A peak list table with a clustered column

- NMRExperiments:

  Two and only two experiments to compare in the plot

- chemshift_range:

  A region, make it so it does not cover a huge range (maybe 1ppm or
  less)

- baselineThresh:

  If given (as returned from the
  [`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md))
  the baseline threshold will be plotted. This can be useful to diagnose
  whether a peak is missing due to this threshold or due to other
  parameters (e.g. `SNR.Th`). See
  [`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md).

## Value

A plot of the two experiments in the given chemshift range, with lines
connecting peaks identified as the same and dots showing peaks without
pairs
