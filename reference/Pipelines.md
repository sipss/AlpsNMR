# Pipelines

Uses
[nmr_pca_outliers_robust](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_robust.md)
to perform the detection of outliers

Normalize the full spectra to the internal calibrant region, then
exclude that region and finally perform PQN normalization.

## Usage

``` r
pipe_load_samples(samples_dir, glob = "*0", output_dir = NULL)

pipe_add_metadata(nmr_dataset_rds, excel_file, output_dir)

pipe_interpolate_1D(nmr_dataset_rds, axis, output_dir)

pipe_exclude_regions(nmr_dataset_rds, exclude, output_dir)

pipe_outlier_detection(nmr_dataset_rds, output_dir)

pipe_filter_samples(nmr_dataset_rds, conditions, output_dir)

pipe_peakdet_align(
  nmr_dataset_rds,
  nDivRange_ppm = 0.1,
  scales = seq(1, 16, 2),
  baselineThresh = 0.01,
  SNR.Th = -1,
  maxShift_ppm = 0.0015,
  acceptLostPeak = FALSE,
  output_dir = NULL
)

pipe_peak_integration(
  nmr_dataset_rds,
  peak_det_align_dir,
  peak_width_ppm,
  output_dir
)

pipe_normalization(
  nmr_dataset_rds,
  internal_calibrant = NULL,
  output_dir = NULL
)
```

## Arguments

- samples_dir:

  The directory where the samples are

- glob:

  A wildcard aka globbing pattern (e.g. `*.csv`) passed on to
  [`grep()`](https://rdrr.io/r/base/grep.html) to filter paths.

- output_dir:

  The output directory for this pipe element

- nmr_dataset_rds:

  The nmr_dataset.rds file name coming from previous nodes

- excel_file:

  An excel file name. See details for the requirements

  The excel file can have one or more sheets. The excel sheets need to
  be as simple as possible: One header column on the first row and
  values below.

  Each of the sheets contain metadata that has to be integrated. The
  merge (technically a left join) is done using the first column of each
  sheet as key.

  In practical terms this means that the first sheet of the excel file
  MUST start with an "NMRExperiment" column, and as many additional
  columns to add (e.g. FluidXBarcode, SampleCollectionDate, TimePoint
  and SubjectID).

  The second sheet can have as the first column any of the already added
  columns, for instance the "SubjectID", and any additional columns
  (e.g. Gender, Age).

  The first column on each sheet, named the key column, MUST have unique
  values. For instance, a sheet starting with "SubjectID" MUST specify
  each subject ID only once (without repetitions).

- axis:

  The ppm axis range and optionally the ppm step. Set it to `NULL` for
  autodetection

- exclude:

  A list with regions to be removed Typically:
  `exclude = list(water = c(4.7, 5.0))`

- conditions:

  A character vector with conditions to filter metadata. The
  `conditions` parameter should be a character vector of valid R logical
  conditions. Some examples:

  - conditions \<- 'Gender == "Female"'

  - conditions \<- 'Cohort == "Chuv"'

  - conditions \<- 'TimePoint %in% c("T0", "T31")'

  - conditions \<- c(Cohort == "Chuv", 'TimePoint %in% c("T0", "T31")')

  Only samples fullfilling all the given conditions are kept in further
  analysis.

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
    ([`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md))
    with the given `range_without_peaks`. There is no ppm range
    guaranteed to be free of peaks for every sample type, so
    `range_without_peaks` must be given when `baselineThresh` is `NULL`.

- SNR.Th:

  The parameter of peakDetectionCWT function of MassSpecWavelet package,
  look it up in the original function. If you set -1, the function will
  itself re-compute this value.

- maxShift_ppm:

  The maximum shift allowed, in ppm

- acceptLostPeak:

  This is an option for users, TRUE is the default value. If the users
  believe that all the peaks in the peak list are true positive, change
  it to FALSE.

- peak_det_align_dir:

  Output directory from pipe_peakdet_align

- peak_width_ppm:

  A peak width in ppm

- internal_calibrant:

  A ppm range where the internal calibrant is, or `NULL`.

## Value

This function saves the result to the output directory

This function saves the result to the output directory

This function saves the result to the output directory

This function saves the result to the output directory

This function saves the result to the output directory

Pipeline: Filter samples according to metadata conditions

Pipeline: Peak detection and Alignment

Pipeline: Peak integration

Pipe: Full spectra normalization

## Details

If there is no internal calibrant, only the PQN normalization is done.

## See also

Other import/export functions:
[`files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/files_to_rDolphin.md),
[`load_and_save_functions`](https://sipss.github.io/AlpsNMR/reference/load_and_save_functions.md),
[`nmr_data()`](https://sipss.github.io/AlpsNMR/reference/nmr_data.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_read_bruker_fid()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_bruker_fid.md),
[`nmr_read_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_samples.md),
[`nmr_zip_bruker_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_zip_bruker_samples.md),
[`save_files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/save_files_to_rDolphin.md),
[`save_profiling_output()`](https://sipss.github.io/AlpsNMR/reference/save_profiling_output.md),
[`to_ChemoSpec()`](https://sipss.github.io/AlpsNMR/reference/to_ChemoSpec.md)

Other metadata functions:
[`nmr_meta_add()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_add.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_meta_get()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get.md),
[`nmr_meta_get_column()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get_column.md),
[`nmr_meta_groups()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_groups.md)

Other outlier detection functions:
[`nmr_pca_outliers()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers.md),
[`nmr_pca_outliers_filter()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_filter.md),
[`nmr_pca_outliers_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_plot.md),
[`nmr_pca_outliers_robust()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_robust.md)

Other peak detection functions:
[`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md),
[`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md),
[`nmr_detect_peaks_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot.md),
[`nmr_detect_peaks_plot_overview()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot_overview.md),
[`nmr_detect_peaks_tune_snr()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_tune_snr.md),
[`nmr_identify_regions_blood()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_blood.md),
[`nmr_identify_regions_cell()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_cell.md),
[`nmr_identify_regions_urine()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_urine.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)

Other alignment functions:
[`nmr_align()`](https://sipss.github.io/AlpsNMR/reference/nmr_align.md),
[`nmr_align_find_ref()`](https://sipss.github.io/AlpsNMR/reference/nmr_align_find_ref.md)

Other peak integration functions:
[`get_integration_with_metadata()`](https://sipss.github.io/AlpsNMR/reference/get_integration_with_metadata.md),
[`nmr_identify_regions_blood()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_blood.md),
[`nmr_identify_regions_cell()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_cell.md),
[`nmr_identify_regions_urine()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_urine.md),
[`nmr_integrate_peak_positions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_peak_positions.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)

## Examples

``` r
## Example of pipeline usage
## There are differet ways of load the dataset
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
# excel_file <- system.file("dataset-demo",
#                          "dummy_metadata.xlsx",
#                          package = "AlpsNMR")
# output_dir <- tempdir()

## Load samples with pipes
# pipe_load_samples(dir_to_demo_dataset,
#                  glob = "*.zip",
#                  output_dir = "../pipe_output")

## Another way to load it
# nmr_dataset <- nmr_read_samples_dir(dir_to_demo_dataset)

## Saving the dataset in a .rds file
# nmr_dataset_rds <- tempfile(fileext = ".rds")
# nmr_dataset_save(nmr_dataset, nmr_dataset_rds)

## Interpolation
# pipe_interpolate_1D(nmr_dataset_rds,
#                    axis = c(min = -0.5, max = 10, by = 2.3E-4),
#                    output_dir)

## Get the new path, based in output_dir
# nmr_dataset_rds <- paste(output_dir, "\", "nmr_dataset.rds", sep = "", collapse = NULL)

## Adding metadata to samples
# pipe_add_metadata(nmr_dataset_rds = nmr_dataset_rds, output_dir = output_dir,
#                  excel_file = excel_file)

## Filtering samples
# conditions <- 'SubjectID == "Ana"'
# pipe_filter_samples(nmr_dataset_rds, conditions, output_dir)

## Outlier detection
# pipe_outlier_detection(nmr_dataset_rds, output_dir)

## Exclude regions
# exclude_regions <- list(water = c(5.1, 4.5))
# pipe_exclude_regions(nmr_dataset_rds, exclude_regions, output_dir)

## peak aling
# pipe_peakdet_align(nmr_dataset_rds, output_dir = output_dir)

## peak integration
# pipe_peak_integration(nmr_dataset_rds,
#                      peak_det_align_dir = output_dir,
#                      peak_width_ppm = 0.006, output_dir)

## Normalization
# pipe_normalization(nmr_dataset_rds, output_dir = output_dir)
```
