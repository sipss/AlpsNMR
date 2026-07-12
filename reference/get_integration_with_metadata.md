# Get integrals with metadata from `integrate peak positions`

Get integrals with metadata from `integrate peak positions`

## Usage

``` r
get_integration_with_metadata(integration_object)
```

## Arguments

- integration_object:

  A
  [nmr_dataset](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md)
  object

## Value

Get integrals with metadata from `integrate peak positions`

integration dataframe

## See also

Other peak integration functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`nmr_identify_regions_blood()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_blood.md),
[`nmr_identify_regions_cell()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_cell.md),
[`nmr_identify_regions_urine()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_urine.md),
[`nmr_integrate_peak_positions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_peak_positions.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)

Other nmr_dataset_1D functions:
[`[.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset_1D.md),
[`format.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_1D.md),
[`is.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_1D.md),
[`nmr_integrate_peak_positions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_peak_positions.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md),
[`nmr_meta_add()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_add.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_meta_get()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get.md),
[`nmr_meta_get_column()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get_column.md),
[`nmr_ppm_resolution()`](https://sipss.github.io/AlpsNMR/reference/nmr_ppm_resolution.md),
[`print.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_1D.md)

## Examples

``` r
peak_table <- matrix(1:6, nrow = 2, ncol = 3)
rownames(peak_table) <- c("10", "20")
colnames(peak_table) <- c("ppm_1.2", "ppm1.4", "ppm1.6")

dataset <- new_nmr_dataset_peak_table(
    peak_table = peak_table,
    metadata = list(external = data.frame(NMRExperiment = c("10", "20"), Condition = c("A", "B")))
)
get_integration_with_metadata(dataset)
#> # A tibble: 2 × 5
#>   NMRExperiment Condition ppm_1.2 ppm1.4 ppm1.6
#>   <chr>         <chr>       <int>  <int>  <int>
#> 1 10            A               1      3      5
#> 2 20            B               2      4      6
```
