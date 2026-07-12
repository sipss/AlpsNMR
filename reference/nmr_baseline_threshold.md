# Threshold estimation for peak detection

Estimates the threshold value for peak detection on an
[nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
object by examining a range without peaks that you must provide (there
is no ppm range guaranteed to be free of peaks for every sample type).

## Usage

``` r
nmr_baseline_threshold(
  nmr_dataset,
  range_without_peaks = NULL,
  method = c("mean3sd", "median3mad")
)
```

## Arguments

- nmr_dataset:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md).

- range_without_peaks:

  A vector with two doubles describing a range without peaks suitable
  for baseline detection. There is no such a range that works for every
  sample type, so you must inspect your spectra and provide one.

- method:

  Either "mean3sd" or the more robust "median3mad". See the details.

## Value

Numerical. A threshold value in intensity below that no peak is
detected.

## Details

Two methods can be used:

- "mean3sd": The mean3sd method computes the mean and the standard
  deviation of each spectrum in the given range. The mean spectrum and
  the mean standard deviation are both vectors of length equal to the
  number of points in the given range. The mean of the mean spectrum the
  noise. The threshold is defined as `center + 3 dispersion`, and it is
  one single threshold for the whole dataset. This is the default for
  backwards compatibility.

- "median3mad": First we take the data matrix. If we have estimated a
  baseline already, subtract it. In the defined region without peaks,
  estimate the median of each sample and its median absolute deviation.
  Return a vector of length equal to the number of samples with the
  `median+3mad` for each sample. This is a new more robust method.

## See also

Other peak detection functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md),
[`nmr_detect_peaks_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot.md),
[`nmr_detect_peaks_plot_overview()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot_overview.md),
[`nmr_detect_peaks_tune_snr()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_tune_snr.md),
[`nmr_identify_regions_blood()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_blood.md),
[`nmr_identify_regions_cell()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_cell.md),
[`nmr_identify_regions_urine()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_urine.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)

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
```
