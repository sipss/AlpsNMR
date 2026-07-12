# Integrate peak positions

The function allows the integration of a given ppm vector with a
specific width.

## Usage

``` r
nmr_integrate_peak_positions(
  samples,
  peak_pos_ppm,
  peak_width_ppm = 0.006,
  ...
)
```

## Arguments

- samples:

  A
  [nmr_dataset](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md)
  object

- peak_pos_ppm:

  The peak positions, in ppm

- peak_width_ppm:

  The peak widths (or a single peak width for all peaks)

- ...:

  Arguments passed on to
  [`nmr_integrate_regions`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)

  `regions`

  : A named list. Each element of the list is a region, given as a named
    numeric vector of length two with the range to integrate. The name
    of the region will be the name of the column

  `fix_baseline`

  : A logical. If `TRUE` it removes the baseline. See details below

  `excluded_regions_as_zero`

  : A logical. It determines the behaviour of the integration when
    integrating regions that have been excluded. If `TRUE`, it will
    treat those regions as zero. If `FALSE` (the default) it will return
    NA values.

    If `fix_baseline` is `TRUE`, then the region boundaries are used to
    estimate a baseline. The baseline is estimated "connecting the
    boundaries with a straight line". Only when the spectrum is above
    the baseline the area is integrated (negative contributions due to
    the baseline estimation are ignored).

  `set_negative_areas_to_zero`

  : A logical. Ignored if `fix_baseline` is `FALSE`. When set to `TRUE`
    negative areas are set to zero.

## Value

Integrate peak positions

## See also

Other peak integration functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`get_integration_with_metadata()`](https://sipss.github.io/AlpsNMR/reference/get_integration_with_metadata.md),
[`nmr_identify_regions_blood()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_blood.md),
[`nmr_identify_regions_cell()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_cell.md),
[`nmr_identify_regions_urine()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_urine.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)

Other nmr_dataset_1D functions:
[`[.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset_1D.md),
[`format.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_1D.md),
[`get_integration_with_metadata()`](https://sipss.github.io/AlpsNMR/reference/get_integration_with_metadata.md),
[`is.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_1D.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md),
[`nmr_meta_add()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_add.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_meta_get()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get.md),
[`nmr_meta_get_column()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get_column.md),
[`nmr_ppm_resolution()`](https://sipss.github.io/AlpsNMR/reference/nmr_ppm_resolution.md),
[`print.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_1D.md)
