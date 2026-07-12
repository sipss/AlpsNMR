# Export 1D NMR data to a CSV file

Export 1D NMR data to a CSV file

## Usage

``` r
nmr_export_data_1r(nmr_dataset, filename)
```

## Arguments

- nmr_dataset:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object

- filename:

  The csv filename

## Value

The nmr_dataset object (unmodified)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
# nmr_export_data_1r(dataset_1D, "exported_nmr_dataset")
```
