# Print for nmr_dataset

Print for nmr_dataset

## Usage

``` r
# S3 method for class 'nmr_dataset'
print(x, ...)
```

## Arguments

- x:

  an
  [nmr_dataset](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md)
  object

- ...:

  for future use

## Value

Print for nmr_dataset

## See also

Other class helper functions:
[`format.nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset.md),
[`format.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_1D.md),
[`format.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_peak_table.md),
[`is.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_1D.md),
[`is.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_peak_table.md),
[`new_nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/new_nmr_dataset.md),
[`new_nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/new_nmr_dataset_1D.md),
[`new_nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/new_nmr_dataset_peak_table.md),
[`print.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_1D.md),
[`print.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_peak_table.md),
[`validate_nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset.md),
[`validate_nmr_dataset_family()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset_family.md),
[`validate_nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset_peak_table.md)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
print(dataset)
#> An nmr_dataset (3 samples) 
```
