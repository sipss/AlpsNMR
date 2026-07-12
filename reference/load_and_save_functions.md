# nmr_dataset_load

nmr_dataset_load

nmr_dataset_save

## Usage

``` r
nmr_dataset_load(file_name)

nmr_dataset_save(nmr_dataset, file_name, ...)
```

## Arguments

- file_name:

  The file name to load or save to

- nmr_dataset:

  An object from the
  [nmr_dataset_family](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_family.md)

- ...:

  Additional arguments passed to
  [saveRDS](https://rdrr.io/r/base/readRDS.html).

## Value

Functions to load and save nmr_dataset objects

load nmr dataset

save nmr dataset

## See also

Other import/export functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/files_to_rDolphin.md),
[`nmr_data()`](https://sipss.github.io/AlpsNMR/reference/nmr_data.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_read_bruker_fid()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_bruker_fid.md),
[`nmr_read_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_samples.md),
[`nmr_zip_bruker_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_zip_bruker_samples.md),
[`save_files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/save_files_to_rDolphin.md),
[`save_profiling_output()`](https://sipss.github.io/AlpsNMR/reference/save_profiling_output.md),
[`to_ChemoSpec()`](https://sipss.github.io/AlpsNMR/reference/to_ChemoSpec.md)

## Examples

``` r
# dataset <- nmr_dataset_load("test")
nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))

dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
# nmr_dataset_save(dataset, "test")
```
