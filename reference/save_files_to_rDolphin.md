# Save files to rDoplhin

The function saves the CSV files required by to_rDolphin and
Automatic_targeted_profiling functions for metabolite profiling.

## Usage

``` r
save_files_to_rDolphin(files_rDolphin, output_directory)
```

## Arguments

- files_rDolphin:

  a list containing 4 elements from `files_to_rDolphin`

  - `meta_rDolphin`: metadata in rDolphin format,

  - `NMR_spectra`: spectra matrix

  - `ROI`: ROI template

  - `Parameters_blood`: parameters file

- output_directory:

  a directory in which the CSV files are saved

## Value

CSV files containing:

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
[`save_profiling_output()`](https://sipss.github.io/AlpsNMR/reference/save_profiling_output.md),
[`to_ChemoSpec()`](https://sipss.github.io/AlpsNMR/reference/to_ChemoSpec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dataset <- system.file("dataset-demo", package = "AlpsNMR")
excel_file <- system.file("dataset-demo", "dummy_metadata.xlsx", package = "AlpsNMR")
nmr_dataset <- nmr_read_samples_dir(dataset)
files_rDolphin <- files_to_rDolphin_blood(nmr_dataset)
save_files_to_rDolphin(files_rDolphin, output_directory = ".")
} # }
```
