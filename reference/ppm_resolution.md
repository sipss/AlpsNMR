# Unlisted PPM resolution

A wrapper to unlist the output from the function
`nmr_ppm_resolution(nmr_dataset)` when no interpolation has been
applied.

## Usage

``` r
ppm_resolution(nmr_dataset)
```

## Arguments

- nmr_dataset:

  An object containing NMR samples

## Value

A number (the ppm resolution, measured in ppms)

Numeric (the ppm resolution, measured in ppms)

## Examples

``` r
nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
nmr_ppm_resolution(nmr_dataset)
#> [1] 0.02
```
