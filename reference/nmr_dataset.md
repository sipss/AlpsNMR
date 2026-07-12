# nmr_dataset (S3 class)

An `nmr_dataset` represents a set of NMR samples. It is defined as an S3
class, and it can be treated as a regular list.

## Details

It currently has the following elements:

- `metadata`: A list of data frames. Each data frame contains metadata
  of a given area (acquisition parameters, preprocessing parameters,
  general sample information...)

- `axis`: A list with length equal to the dimensionality of the data.
  For 1D spectra it is a list with a numeric vector

- `data_*`: Data arrays with the actual spectra. The first index
  represents the sample, the rest of the indices match the length of
  each `axis`. Typically `data_1r` is a matrix with one sample on each
  row and the chemical shifts in the columns.

- `num_samples`: The number of samples in the dataset

## See also

[Functions to save and load these
objects](https://sipss.github.io/AlpsNMR/reference/load_and_save_functions.md)

Other AlpsNMR dataset objects:
[`nmr_dataset_family`](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_family.md)

## Examples

``` r
metadata_1D <- list(external = data.frame(NMRExperiment = c("10", "20")))
# Sample 10 and Sample 20 can have different lengths (due to different setups)
data_fields_1D <- list(data_1r = list(runif(16), runif(32)))
# Each sample has its own axis list, with one element (because this example is 1D)
axis_1D <- list(list(1:16), list(1:32))
my_1D_data <- new_nmr_dataset(metadata_1D, data_fields_1D, axis_1D)
```
