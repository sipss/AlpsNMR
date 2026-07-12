# Plot peak detection results

Plot peak detection results

## Usage

``` r
nmr_detect_peaks_plot(
  nmr_dataset,
  peak_data,
  NMRExperiment = NULL,
  peak_id = NULL,
  accepted_only = NULL,
  ...
)
```

## Arguments

- nmr_dataset:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md).

- peak_data:

  The peak table returned by
  [nmr_detect_peaks](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md)

- NMRExperiment:

  a single NMR experiment to plot

- peak_id:

  A character vector. If given, plot only that peak id.

- accepted_only:

  If `peak_data` contains a logical column named `accepted`, only those
  with `accepted=TRUE` will be counted. By default,
  `accepted_only = TRUE`, unless a `peak_id` is given

- ...:

  Arguments passed to
  [plot.nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/plot.nmr_dataset_1D.md)
  (`chemshift_range`, `...`)

## Value

Plot peak detection results

## See also

Peak_detection nmr_detect_peaks

Other peak detection functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md),
[`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md),
[`nmr_detect_peaks_plot_overview()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot_overview.md),
[`nmr_detect_peaks_tune_snr()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_tune_snr.md),
[`nmr_identify_regions_blood()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_blood.md),
[`nmr_identify_regions_cell()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_cell.md),
[`nmr_identify_regions_urine()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_urine.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)
