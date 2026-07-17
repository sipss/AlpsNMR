# Diagnose SNR threshold in peak detection

Diagnose SNR threshold in peak detection

## Usage

``` r
nmr_detect_peaks_tune_snr(
  ds,
  NMRExperiment = NULL,
  SNR_thresholds = seq(from = 2, to = 6, by = 0.1),
  ...
)
```

## Arguments

- ds:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  dataset

- NMRExperiment:

  A string with the single NMRExperiment used explore the SNR
  thresholds. If not given, use the first one.

- SNR_thresholds:

  A numeric vector with the SNR thresholds to explore

- ...:

  Arguments passed on to
  [`nmr_detect_peaks`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md)

  `nmr_dataset`

  : An
    [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md).

  `nDivRange_ppm`

  : Segment size, in ppms, to divide the spectra and search for peaks.

  `baselineThresh`

  : All peaks with intensities below the thresholds are excluded.
    Either:

    - A numeric vector of length the number of samples. Each number is a
      threshold for that sample

    - A single number. All samples use this number as baseline
      threshold.

    - `NULL`. If that's the case, a default function is used
      ([`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md))
      with the given `range_without_peaks`. There is no ppm range
      guaranteed to be free of peaks for every sample type, so
      `range_without_peaks` must be given when `baselineThresh` is
      `NULL`.

  `range_without_peaks`

  : A numeric vector of length two with a region without peaks. Required
    when `baselineThresh = NULL`, ignored otherwise.

  `fit_lorentzians`

  : If `TRUE`, fit a lorentzian to each detected peak, to infer its
    inflection points. For now disabled for backwards compatibility.

  `verbose`

  : Logical (`TRUE` or `FALSE`). Show informational messages, such as
    the estimated baseline

  `scales`

  : The parameter of peakDetectionCWT function of MassSpecWavelet
    package, look it up in the original function.

  `SNR.Th`

  : The parameter of peakDetectionCWT function of MassSpecWavelet
    package, look it up in the original function. If you set -1, the
    function will itself re-compute this value.

## Value

A list with the following elements:

- `peaks_detected`: A data frame with the columns from the
  [nmr_detect_peaks](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md)
  output and an additional column `SNR_threshold` with the threshold
  used on each row.

- `num_peaks_per_region`: A summary of the `peaks_detected` table, with
  the number of peaks detected on each chemical shift region

- `plot_num_peaks_per_region`: A visual representation of
  `num_peaks_per_region`

- `plot_spectrum_and_detections`: A visual representation of the
  spectrum and the peaks detected with each SNR threshold. Use
  [plotly::ggplotly](https://rdrr.io/pkg/plotly/man/ggplotly.html) or
  [plot_interactive](https://sipss.github.io/AlpsNMR/reference/plot_interactive.md)
  on this to zoom and explore the results.

## See also

nmr_detect_peaks

Other peak detection functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md),
[`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md),
[`nmr_detect_peaks_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot.md),
[`nmr_detect_peaks_plot_overview()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot_overview.md),
[`nmr_identify_regions_blood()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_blood.md),
[`nmr_identify_regions_cell()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_cell.md),
[`nmr_identify_regions_urine()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_urine.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)
