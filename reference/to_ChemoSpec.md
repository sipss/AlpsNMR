# Convert to ChemoSpec Spectra class

Convert to ChemoSpec Spectra class

## Usage

``` r
to_ChemoSpec(nmr_dataset, desc = "A nmr_dataset", group = NULL)
```

## Arguments

- nmr_dataset:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object

- desc:

  a description for the dataset

- group:

  A string with the column name from the metadata that has grouping
  information

## Value

A Spectra object from the ChemoSpec package

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
[`save_profiling_output()`](https://sipss.github.io/AlpsNMR/reference/save_profiling_output.md)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
chemo_spectra <- to_ChemoSpec(dataset_1D)
```
