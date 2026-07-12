# Save rDoplhin output

The function saves the output from Automatic_targeted_profiling function
in CSV format.

## Usage

``` r
save_profiling_output(targeted_profiling, output_directory)
```

## Arguments

- targeted_profiling:

  A list from Automatic_targeted_profiling function

- output_directory:

  a directory in which the CSV files are saved

## Value

rDolphin output from Automatic_targeted_profiling function:

- metabolites_intensity

- metabolites_quantification

- ROI_profiles_used

- chemical_shift

- fitting_error

- half_bandwidth

- signal_area_ratio

## See also

Other import/export functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/files_to_rDolphin.md),
[`load_and_save_functions`](https://sipss.github.io/AlpsNMR/reference/load_and_save_functions.md),
[`nmr_data()`](https://sipss.github.io/AlpsNMR/reference/nmr_data.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_read_bruker_fid()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_bruker_fid.md),
[`nmr_read_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_samples.md),
[`nmr_zip_bruker_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_zip_bruker_samples.md),
[`save_files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/save_files_to_rDolphin.md),
[`to_ChemoSpec()`](https://sipss.github.io/AlpsNMR/reference/to_ChemoSpec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
rDolphin_object <- to_rDolphin(parameters)
targeted_profiling <- Automatic_targeted_profiling(rDolphin)
save_profiling_output(targeted_profiling, output_directory)
} # }
```
