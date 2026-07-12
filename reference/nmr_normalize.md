# Normalize nmr_dataset_1D samples

The `nmr_normalize` function is used to normalize all the samples
according to a given criteria.

## Usage

``` r
nmr_normalize(
  samples,
  method = c("area", "max", "value", "region", "pqn", "none"),
  ...
)

nmr_normalize_extra_info(samples)
```

## Arguments

- samples:

  A
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object

- method:

  The criteria to be used for normalization

  - area: Normalize to the total area

  - max: Normalize to the maximum intensity

  - value: Normalize each sample to a user defined value

  - region: Integrate a region and normalize each sample to that region

  - pqn: Use Probabalistic Quotient Normalization for normalization

  - none: Do not normalize at all

- ...:

  Method dependent arguments:

  - `method == "value"`: - `value`: A numeric vector with the
    normalization values to use

  - `method == "region"`: - `ppm_range`: A chemical shift region to
    integrate - `...`: Other arguments passed on to
    [nmr_integrate_regions](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)

## Value

The
[nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
object, with the samples normalized. Further information for diagnostic
of the normalization process is also saved and can be extracted by
calling `nmr_normalize_extra_info()` afterwards.

## Details

The aim is to correct from changes between samples, so no matter the
criteria used to normalize, once we get the factors (e.g. the areas), we
divide them by the median normalization factor to avoid introducing
global scaling factors.

The `nmr_normalize_extra_info` function is used to extract additional
information after the normalization. Typically, we want to know what was
the actual normalization factor applied to each sample. The extra
information includes a plot, representing the dispersion of the
normalization factor for each sample.

## See also

Other basic functions:
[`nmr_exclude_region()`](https://sipss.github.io/AlpsNMR/reference/nmr_exclude_region.md)

## Examples

``` r
nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
nmr_dataset <- nmr_normalize(nmr_dataset, method = "area")
norm_dataset <- nmr_normalize(nmr_dataset)
norm_dataset$plot
#> NULL
nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
nmr_dataset <- nmr_normalize(nmr_dataset, method = "area")
norm_extra_info <- nmr_normalize_extra_info(nmr_dataset)
norm_extra_info$plot
```
