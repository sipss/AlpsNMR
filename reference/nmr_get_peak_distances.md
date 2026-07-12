# Compute peak to peak distances

Compute peak to peak distances

## Usage

``` r
nmr_get_peak_distances(peak_data, same_sample_dist_factor = 3)
```

## Arguments

- peak_data:

  A peak list

- same_sample_dist_factor:

  The distance between two peaks from the same sample are set to this
  factor multiplied by the maximum of all the peak distances

## Value

A dist object with the peak2peak distances

## Examples

``` r
peak_data <- data.frame(
    NMRExperiment = c("10", "10", "20", "20"),
    peak_id = paste0("Peak", 1:4),
    ppm = c(1, 2, 1.1, 3)
)
peak2peak_dist <- nmr_get_peak_distances(peak_data)
stopifnot(abs(as.numeric(peak2peak_dist) - c(6, 0.1, 2, 0.9, 1, 6)) < 1E-8)
```
