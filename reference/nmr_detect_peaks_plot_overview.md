# Overview of the peak detection results

This plot allows to explore the performance of the peak detection across
all the samples, by summarizing how many peaks are detected on each
sample at each chemical shift range.

## Usage

``` r
nmr_detect_peaks_plot_overview(
  peak_data,
  ppm_breaks = pretty(range(peak_data$ppm), n = 20),
  accepted_only = TRUE
)
```

## Arguments

- peak_data:

  The output of
  [`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md)

- ppm_breaks:

  A numeric vector with the breaks that will be used to count the number
  of the detected peaks.

- accepted_only:

  If `peak_data` contains a logical column named `accepted`, only those
  with `accepted=TRUE` will be counted.

## Value

A scatter plot, with samples on one axis and chemical shift bins in the
other axis. The size of each dot represents the number of peaks found on
a sample within a chemical shift range.

## Details

You can use this plot to find differences in the number of detected
peaks across your dataset, and then use
[`nmr_detect_peaks_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot.md)
to have a finer look at specific samples and chemical shifts, and assess
graphically that the peak detection results that you have are correct.

## See also

Peak_detection

Other peak detection functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md),
[`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md),
[`nmr_detect_peaks_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot.md),
[`nmr_detect_peaks_tune_snr()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_tune_snr.md),
[`nmr_identify_regions_blood()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_blood.md),
[`nmr_identify_regions_cell()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_cell.md),
[`nmr_identify_regions_urine()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_urine.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)
