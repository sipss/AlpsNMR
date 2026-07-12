# Exclude region from samples

Excludes a given region (for instance to remove the water peak)

## Usage

``` r
nmr_exclude_region(samples, exclude = list(water = c(4.7, 5)))

# S3 method for class 'nmr_dataset_1D'
nmr_exclude_region(samples, exclude = list(water = c(4.7, 5)))
```

## Arguments

- samples:

  An object

- exclude:

  A list with regions to be removed Typically:
  `exclude = list(water = c(4.7, 5.0))`

## Value

The same object, with the regions excluded

## See also

Other basic functions:
[`nmr_normalize()`](https://sipss.github.io/AlpsNMR/reference/nmr_normalize.md)

## Examples

``` r
nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
exclude_regions <- list(water = c(5.1, 4.5))
nmr_dataset <- nmr_exclude_region(nmr_dataset, exclude = exclude_regions)

nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
exclude_regions <- list(water = c(5.1, 4.5))
nmr_dataset <- nmr_exclude_region(nmr_dataset, exclude = exclude_regions)
```
