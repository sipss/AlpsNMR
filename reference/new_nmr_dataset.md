# Create an nmr_dataset object

Create an nmr_dataset object

## Usage

``` r
new_nmr_dataset(metadata, data_fields, axis)
```

## Arguments

- metadata:

  A named list of data frames

- data_fields:

  A named list. Check the examples

- axis:

  A list. Check the examples

## Value

Create an nmr_dataset object

Create an nmr_dataset object

## See also

Other class helper functions:
[`format.nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset.md),
[`format.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_1D.md),
[`format.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_peak_table.md),
[`is.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_1D.md),
[`is.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_peak_table.md),
[`new_nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/new_nmr_dataset_1D.md),
[`new_nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/new_nmr_dataset_peak_table.md),
[`print.nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset.md),
[`print.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_1D.md),
[`print.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_peak_table.md),
[`validate_nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset.md),
[`validate_nmr_dataset_family()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset_family.md),
[`validate_nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset_peak_table.md)

## Examples

``` r
#
metadata_1D <- list(external = data.frame(NMRExperiment = c("10", "20")))
# Sample 10 and Sample 20 can have different lengths (due to different setups)
data_fields_1D <- list(data_1r = list(runif(16), runif(32)))
# Each sample has its own axis list, with one element (because this example is 1D)
axis_1D <- list(list(1:16), list(1:32))
my_1D_data <- new_nmr_dataset(metadata_1D, data_fields_1D, axis_1D)

# Example for 2D samples
metadata_2D <- list(external = data.frame(NMRExperiment = c("11", "21")))
data_fields_2D <- list(data_2rr = list(matrix(runif(16 * 3), nrow = 16, ncol = 3),
    runif(32 * 3),
    nrow = 32, ncol = 3
))
# Each sample has its own axis list, with one element (because this example is 1D)
axis_2D <- list(list(1:16, 1:3), list(1:32, 1:3))
my_2D_data <- new_nmr_dataset(metadata_2D, data_fields_2D, axis_2D)
```
