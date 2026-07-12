# Align NMR spectra

This function is based on
[speaq::dohCluster](https://rdrr.io/pkg/speaq/man/dohCluster.html).

## Usage

``` r
nmr_align(
  nmr_dataset,
  peak_data,
  NMRExp_ref = NULL,
  maxShift_ppm = 0.0015,
  acceptLostPeak = FALSE
)
```

## Arguments

- nmr_dataset:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)

- peak_data:

  The detected peak data given by
  [nmr_detect_peaks](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md).

- NMRExp_ref:

  NMRExperiment of the reference to use for alignment

- maxShift_ppm:

  The maximum shift allowed, in ppm

- acceptLostPeak:

  This is an option for users, TRUE is the default value. If the users
  believe that all the peaks in the peak list are true positive, change
  it to FALSE.

## Value

An
[nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md),
with the spectra aligned

## See also

Other alignment functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`nmr_align_find_ref()`](https://sipss.github.io/AlpsNMR/reference/nmr_align_find_ref.md)

Other peak alignment functions:
[`nmr_align_find_ref()`](https://sipss.github.io/AlpsNMR/reference/nmr_align_find_ref.md)
