# Interpolate a set of 1D NMR Spectra

Interpolate a set of 1D NMR Spectra

## Usage

``` r
nmr_interpolate_1D(samples, axis = c(min = 0.2, max = 10, by = 8e-04))

# S3 method for class 'nmr_dataset'
nmr_interpolate_1D(samples, axis = c(min = 0.2, max = 10, by = 8e-04))
```

## Arguments

- samples:

  An NMR dataset

- axis:

  The ppm axis range and optionally the ppm step. Set it to `NULL` for
  autodetection

## Value

Interpolate a set of 1D NMR Spectra

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))

dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
```
