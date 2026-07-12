# Import SummarizedExperiment as 1D NMR data

Import SummarizedExperiment as 1D NMR data

## Usage

``` r
SummarizedExperiment_to_nmr_data_1r(se)
```

## Arguments

- se:

  An SummarizedExperiment object

## Value

nmr_dataset An
[nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
object (unmodified)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
se <- nmr_data_1r_to_SummarizedExperiment(dataset_1D)
dataset_1D <- SummarizedExperiment_to_nmr_data_1r(se)
```
