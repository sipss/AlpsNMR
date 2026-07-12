# Integrate regions

Integrate given regions and return a data frame with them

## Usage

``` r
nmr_integrate_regions(samples, regions, ...)

# S3 method for class 'nmr_dataset_1D'
nmr_integrate_regions(
  samples,
  regions,
  fix_baseline = FALSE,
  excluded_regions_as_zero = FALSE,
  set_negative_areas_to_zero = FALSE,
  ...
)
```

## Arguments

- samples:

  A
  [nmr_dataset](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md)
  object

- regions:

  A named list. Each element of the list is a region, given as a named
  numeric vector of length two with the range to integrate. The name of
  the region will be the name of the column

- ...:

  Keep for compatibility

- fix_baseline:

  A logical. If `TRUE` it removes the baseline. See details below

- excluded_regions_as_zero:

  A logical. It determines the behaviour of the integration when
  integrating regions that have been excluded. If `TRUE`, it will treat
  those regions as zero. If `FALSE` (the default) it will return NA
  values.

  If `fix_baseline` is `TRUE`, then the region boundaries are used to
  estimate a baseline. The baseline is estimated "connecting the
  boundaries with a straight line". Only when the spectrum is above the
  baseline the area is integrated (negative contributions due to the
  baseline estimation are ignored).

- set_negative_areas_to_zero:

  A logical. Ignored if `fix_baseline` is `FALSE`. When set to `TRUE`
  negative areas are set to zero.

## Value

An
[nmr_dataset_peak_table](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_peak_table.md)
object

## See also

Other peak detection functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md),
[`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md),
[`nmr_detect_peaks_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot.md),
[`nmr_detect_peaks_plot_overview()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot_overview.md),
[`nmr_detect_peaks_tune_snr()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_tune_snr.md),
[`nmr_identify_regions_blood()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_blood.md),
[`nmr_identify_regions_cell()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_cell.md),
[`nmr_identify_regions_urine()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_urine.md)

Other peak integration functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`get_integration_with_metadata()`](https://sipss.github.io/AlpsNMR/reference/get_integration_with_metadata.md),
[`nmr_identify_regions_blood()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_blood.md),
[`nmr_identify_regions_cell()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_cell.md),
[`nmr_identify_regions_urine()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_urine.md),
[`nmr_integrate_peak_positions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_peak_positions.md)

Other nmr_dataset_1D functions:
[`[.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset_1D.md),
[`format.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_1D.md),
[`get_integration_with_metadata()`](https://sipss.github.io/AlpsNMR/reference/get_integration_with_metadata.md),
[`is.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_1D.md),
[`nmr_integrate_peak_positions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_peak_positions.md),
[`nmr_meta_add()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_add.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_meta_get()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get.md),
[`nmr_meta_get_column()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get_column.md),
[`nmr_ppm_resolution()`](https://sipss.github.io/AlpsNMR/reference/nmr_ppm_resolution.md),
[`print.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_1D.md)

## Examples

``` r
# Creating a dataset
dataset <- new_nmr_dataset_1D(
    ppm_axis = 1:10,
    data_1r = matrix(sample(0:99, replace = TRUE), nrow = 10),
    metadata = list(external = data.frame(NMRExperiment = c(
        "10",
        "20", "30", "40", "50", "60", "70", "80", "90", "100"
    )))
)

# Integrating selected regions
peak_table_integration <- nmr_integrate_regions(
    samples = dataset,
    regions = list(ppm = c(2, 5))
)

# Creating a dataset
dataset <- new_nmr_dataset_1D(
    ppm_axis = 1:10,
    data_1r = matrix(sample(0:99, replace = TRUE), nrow = 10),
    metadata = list(external = data.frame(NMRExperiment = c(
        "10",
        "20", "30", "40", "50", "60", "70", "80", "90", "100"
    )))
)

# Integrating selected regions
peak_table_integration <- nmr_integrate_regions(
    samples = dataset,
    regions = list(ppm = c(2, 5)),
    fix_baseline = FALSE
)
```
