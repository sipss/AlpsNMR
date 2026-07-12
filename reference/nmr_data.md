# Set/Return the full spectra matrix

Set/Return the full spectra matrix

## Usage

``` r
nmr_data(nmr_dataset, ...)

# S3 method for class 'nmr_dataset_1D'
nmr_data(nmr_dataset, what = "data_1r", ...)

nmr_data(nmr_dataset, ...) <- value

# S3 method for class 'nmr_dataset_1D'
nmr_data(nmr_dataset, what = "data_1r", ...) <- value
```

## Arguments

- nmr_dataset:

  An object from the
  [nmr_dataset_family](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_family.md)
  to get the raw data from

- ...:

  Passed on to methods for compatibility

- what:

  What data do we want to get (default: `data_1r`)

- value:

  A matrix

## Value

a matrix

The given nmr_dataset

## See also

Other import/export functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/files_to_rDolphin.md),
[`load_and_save_functions`](https://sipss.github.io/AlpsNMR/reference/load_and_save_functions.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_read_bruker_fid()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_bruker_fid.md),
[`nmr_read_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_samples.md),
[`nmr_zip_bruker_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_zip_bruker_samples.md),
[`save_files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/save_files_to_rDolphin.md),
[`save_profiling_output()`](https://sipss.github.io/AlpsNMR/reference/save_profiling_output.md),
[`to_ChemoSpec()`](https://sipss.github.io/AlpsNMR/reference/to_ChemoSpec.md)

## Examples

``` r
dataset_rds <- system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR")
dataset_1D <- nmr_dataset_load(dataset_rds)
dataset_data <- nmr_data(dataset_1D)
dataset_rds <- system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR")
dataset_1D <- nmr_dataset_load(dataset_rds)
dataset_1D_data <- nmr_data(dataset_1D)
```
