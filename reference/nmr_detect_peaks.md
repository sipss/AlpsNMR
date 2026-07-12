# Peak detection for NMR

The function detects peaks on an
[nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
object, using
[speaq::detectSpecPeaks](https://rdrr.io/pkg/speaq/man/detectSpecPeaks.html).
`detectSpecPeaks` divides the whole spectra into smaller segments and
uses
[MassSpecWavelet::peakDetectionCWT](https://rdrr.io/pkg/MassSpecWavelet/man/peakDetectionCWT.html)
for peak detection.

## Usage

``` r
nmr_detect_peaks(
  nmr_dataset,
  nDivRange_ppm = 0.1,
  scales = seq(1, 16, 2),
  baselineThresh = NULL,
  SNR.Th = 3,
  range_without_peaks = c(9.5, 10),
  fit_lorentzians = FALSE,
  verbose = FALSE
)
```

## Arguments

- nmr_dataset:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md).

- nDivRange_ppm:

  Segment size, in ppms, to divide the spectra and search for peaks.

- scales:

  The parameter of peakDetectionCWT function of MassSpecWavelet package,
  look it up in the original function.

- baselineThresh:

  All peaks with intensities below the thresholds are excluded. Either:

  - A numeric vector of length the number of samples. Each number is a
    threshold for that sample

  - A single number. All samples use this number as baseline threshold.

  - `NULL`. If that's the case, a default function is used
    ([`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md)),
    which assumes that there is no signal in the region 9.5-10 ppm.

- SNR.Th:

  The parameter of peakDetectionCWT function of MassSpecWavelet package,
  look it up in the original function. If you set -1, the function will
  itself re-compute this value.

- range_without_peaks:

  A numeric vector of length two with a region without peaks, only used
  when `baselineThresh = NULL`

- fit_lorentzians:

  If `TRUE`, fit a lorentzian to each detected peak, to infer its
  inflection points. For now disabled for backwards compatibility.

- verbose:

  Logical (`TRUE` or `FALSE`). Show informational messages, such as the
  estimated baseline

## Value

A data frame with the NMRExperiment, the sample index, the position in
ppm and index and the peak intensity

## Details

Optionally afterwards, the peak apex and the peak inflection points are
used to efficiently adjust a lorentzian to each peak, and compute the
peak area and width, as well as the error of the fit. These peak
features can be used afterwards to reject false detections.

## See also

[nmr_align](https://sipss.github.io/AlpsNMR/reference/nmr_align.md) for
peak alignment with the detected peak table

Peak_detection

Other peak detection functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md),
[`nmr_detect_peaks_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot.md),
[`nmr_detect_peaks_plot_overview()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot_overview.md),
[`nmr_detect_peaks_tune_snr()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_tune_snr.md),
[`nmr_identify_regions_blood()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_blood.md),
[`nmr_identify_regions_cell()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_cell.md),
[`nmr_identify_regions_urine()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_urine.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)
