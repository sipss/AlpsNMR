# Get a tidy data frame from nmr_data object

This dataframe is useful for plotting with ggplot, although it may be
very long and therefore use a lot of RAM.

## Usage

``` r
# S3 method for class 'nmr_dataset_1D'
tidy(
  x,
  NMRExperiment = NULL,
  chemshift_range = NULL,
  columns = character(0L),
  matrix_name = "data_1r",
  axis_name = "axis",
  ...
)
```

## Arguments

- x:

  an
  [`nmr_dataset_1D`](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object

- NMRExperiment:

  A character vector with the NMRExperiments to include. `NULL` means
  all.

- chemshift_range:

  range of the chemical shifts to be included. Can be of length 3 to
  include the resolution in the third element (e.g.
  `c(0.2, 0.8, 0.005)`)

- columns:

  A character vector with the metadata columns to get, use `NULL` to get
  all of them.

- matrix_name:

  A string with the matrix name, typically "data_1r"

- axis_name:

  A string with the axis name, for now "axis" is the only valid option

- ...:

  Ignored

## Value

A data frame with `NMRExperiment`, `chemshift`, `intensity` and any
additional column requested

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -1.0, max = 1.6, by = 2.3E-4))
dummy_metadata <- system.file("dataset-demo", "dummy_metadata.xlsx", package = "AlpsNMR")
NMRExp_SubjID <- readxl::read_excel(dummy_metadata, sheet = 1)
dataset_1D <- nmr_meta_add(dataset_1D, NMRExp_SubjID)
df_for_ggplot <- tidy(dataset_1D, chemshift_range = c(1.2, 1.4), columns = "SubjectID")
head(df_for_ggplot)
#>   NMRExperiment chemshift intensity SubjectID
#> 1            10   1.20018  95475.63       Ana
#> 2            20   1.20018 126803.78       Ana
#> 3            30   1.20018 118041.30      Elia
#> 4            10   1.20041  95387.90       Ana
#> 5            20   1.20041 125392.78       Ana
#> 6            30   1.20041 117378.69      Elia
```
