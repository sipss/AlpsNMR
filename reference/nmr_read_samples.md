# Read NMR samples

These functions load samples from files and return a
[nmr_dataset](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md).

## Usage

``` r
nmr_read_samples_dir(
  samples_dir,
  format = "bruker",
  pulse_sequence = NULL,
  metadata_only = FALSE,
  ...
)

nmr_read_samples(
  sample_names,
  format = "bruker",
  pulse_sequence = NULL,
  metadata_only = FALSE,
  ...
)
```

## Arguments

- samples_dir:

  A directory or directories that contain multiple samples

- format:

  Either "bruker" or "jdx"

- pulse_sequence:

  If it is set to a pulse sequence ("NOESY", "JRES", "CPMG"...) it will
  only load the samples that match that pulse sequence.

- metadata_only:

  A logical, to load only metadata (default: `FALSE`)

- ...:

  Arguments passed on to
  [`read_bruker_pdata`](https://sipss.github.io/AlpsNMR/reference/read_bruker_pdata.md)

  `pdata_file`

  : File name of the binary NMR data to load. Usually "1r". If `NULL`,
    it is autodetected based on the dimension

  `sample_path`

  : A character path of the sample directory

  `pdata_path`

  : Path from `sample_path` to the preprocessed data

  `all_components`

  : If `FALSE` load only the real component. Otherwise load the real and
    imaginary components

  `read_pdata_title`

  : If `TRUE` also reads metadata from pdata title file.

- sample_names:

  A character vector with file or directory names.

## Value

a
[nmr_dataset](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md)
object

## See also

[`read_bruker_pdata()`](https://sipss.github.io/AlpsNMR/reference/read_bruker_pdata.md)

Other import/export functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/files_to_rDolphin.md),
[`load_and_save_functions`](https://sipss.github.io/AlpsNMR/reference/load_and_save_functions.md),
[`nmr_data()`](https://sipss.github.io/AlpsNMR/reference/nmr_data.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_read_bruker_fid()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_bruker_fid.md),
[`nmr_zip_bruker_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_zip_bruker_samples.md),
[`save_files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/save_files_to_rDolphin.md),
[`save_profiling_output()`](https://sipss.github.io/AlpsNMR/reference/save_profiling_output.md),
[`to_ChemoSpec()`](https://sipss.github.io/AlpsNMR/reference/to_ChemoSpec.md)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)

dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
zip_files <- fs::dir_ls(dir_to_demo_dataset, glob = "*.zip")
dataset <- nmr_read_samples(sample_names = zip_files)
```
