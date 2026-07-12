# Read Free Induction Decay file

Reads a Bruker FID file. The sample's `acqus` file is used to determine
how to interpret the raw binary data: `BYTORDA` gives the byte order,
`DTYPA` gives the data type (32-bit integer or 64-bit double), `TD`
gives the number of raw (real+imaginary interleaved) data points
actually acquired, and `SW_h` (the spectral width, in Hz) is used to
build the acquisition time axis.

## Usage

``` r
nmr_read_bruker_fid(sample_name)
```

## Arguments

- sample_name:

  A single sample directory. It must contain an `acqus` file and a `fid`
  file.

## Value

A data frame with columns `time_s` (the acquisition time, in seconds, of
each complex data point) and `fid_complex` (the free induction decay, as
a complex vector). Returns `NULL` if the sample has no `fid` file.

## See also

Other import/export functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/files_to_rDolphin.md),
[`load_and_save_functions`](https://sipss.github.io/AlpsNMR/reference/load_and_save_functions.md),
[`nmr_data()`](https://sipss.github.io/AlpsNMR/reference/nmr_data.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_read_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_samples.md),
[`nmr_zip_bruker_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_zip_bruker_samples.md),
[`save_files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/save_files_to_rDolphin.md),
[`save_profiling_output()`](https://sipss.github.io/AlpsNMR/reference/save_profiling_output.md),
[`to_ChemoSpec()`](https://sipss.github.io/AlpsNMR/reference/to_ChemoSpec.md)

## Examples

``` r
fid <- nmr_read_bruker_fid("sample.fid")
```
