# Create one zip file for each brucker sample path

Create one zip file for each brucker sample path

## Usage

``` r
nmr_zip_bruker_samples(path, workdir, overwrite = FALSE, ...)
```

## Arguments

- path:

  Character vector with sample directories

- workdir:

  Directory to store zip files

- overwrite:

  Should existing zip files be overwritten?

- ...:

  Passed to [utils::zip](https://rdrr.io/r/utils/zip.html)

## Value

A character vector of the same length as path, with the zip file names

## See also

Other import/export functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/files_to_rDolphin.md),
[`load_and_save_functions`](https://sipss.github.io/AlpsNMR/reference/load_and_save_functions.md),
[`nmr_data()`](https://sipss.github.io/AlpsNMR/reference/nmr_data.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_read_bruker_fid()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_bruker_fid.md),
[`nmr_read_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_samples.md),
[`save_files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/save_files_to_rDolphin.md),
[`save_profiling_output()`](https://sipss.github.io/AlpsNMR/reference/save_profiling_output.md),
[`to_ChemoSpec()`](https://sipss.github.io/AlpsNMR/reference/to_ChemoSpec.md)

## Examples

``` r
save_zip_files_to <- tempfile(pattern = "zip_file_storage_")
where_your_samples_are <- tempfile(pattern = "where_your_samples_are")
# prepare sample:
zip::unzip(
  system.file("dataset-demo", "10.zip", package = "AlpsNMR"),
  exdir = where_your_samples_are
)

outpaths <- nmr_zip_bruker_samples(
    list.files(where_your_samples_are, full.names = TRUE),
    workdir = save_zip_files_to
)
```
