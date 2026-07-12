# NMR peak identification (plasma/serum samples)

Identify given regions and return a data frame with plausible
assignations in human plasma/serum samples.

## Usage

``` r
nmr_identify_regions_blood(
  ppm_to_assign,
  num_proposed_compounds = 3,
  verbose = FALSE
)
```

## Arguments

- ppm_to_assign:

  A vector with the ppm regions to assign

- num_proposed_compounds:

  set the number of proposed metabolites sorted by the number times
  reported in the HMDB: `HMDB_blood`.

- verbose:

  Logical value. Set it to TRUE to print additional information

## Value

a data frame with plausible assignations.

## See also

Other peak detection functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md),
[`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md),
[`nmr_detect_peaks_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot.md),
[`nmr_detect_peaks_plot_overview()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot_overview.md),
[`nmr_detect_peaks_tune_snr()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_tune_snr.md),
[`nmr_identify_regions_cell()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_cell.md),
[`nmr_identify_regions_urine()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_urine.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)

Other peak integration functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`get_integration_with_metadata()`](https://sipss.github.io/AlpsNMR/reference/get_integration_with_metadata.md),
[`nmr_identify_regions_cell()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_cell.md),
[`nmr_identify_regions_urine()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_urine.md),
[`nmr_integrate_peak_positions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_peak_positions.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)

## Examples

``` r
# We identify regions from from the corresponding ppm storaged in a vector.
ppm_to_assign <- c(
    4.060960203, 3.048970634, 2.405935596,
    3.24146865, 0.990616851, 1.002075066, 0.955325548
)
identification <- nmr_identify_regions_blood(ppm_to_assign)
```
