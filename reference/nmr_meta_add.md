# Add metadata to an nmr_dataset object

This is useful to add metadata to datasets that can be later used for
plotting spectra or further analysis (PCA...).

## Usage

``` r
nmr_meta_add(nmr_data, metadata, by = "NMRExperiment")

nmr_meta_add_tidy_excel(nmr_data, excel_file)
```

## Arguments

- nmr_data:

  an
  [nmr_dataset_family](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_family.md)
  object

- metadata:

  A data frame with metadata to add

- by:

  A column name of both the `nmr_data$metadata$external` and the
  metadata data.frame. If you want to merge two columns with different
  headers you can use a named character vector
  `c("NMRExperiment" = "ExperimentNMR")` where the left side is the
  column name of the `nmr_data$metadata$external` and the right side is
  the column name of the metadata data frame.

- excel_file:

  Path to a tidy Excel file name. The Excel can consist of multiple
  sheets, that are added sequentially. The first column of the first
  sheet MUST be named as one of the metadata already present in the
  dataset, typically will be "NMRExperiment". The rest of the columns of
  the first sheet can be named at will. Similary, the first column of
  the second sheet must be named as one of the metadata already present
  in the dataset, typically "NMRExperiment" or any of the columns of the
  first sheet. The rest of the columns of the second sheet can be named
  at will. See the package vignette for an example.

## Value

The nmr_dataset_family object with the added metadata

## See also

Other metadata functions:
[`Pipelines`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_meta_get()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get.md),
[`nmr_meta_get_column()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get_column.md),
[`nmr_meta_groups()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_groups.md)

Other nmr_dataset functions:
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_meta_get()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get.md),
[`nmr_meta_get_column()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get_column.md)

Other nmr_dataset_1D functions:
[`[.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset_1D.md),
[`format.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_1D.md),
[`get_integration_with_metadata()`](https://sipss.github.io/AlpsNMR/reference/get_integration_with_metadata.md),
[`is.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_1D.md),
[`nmr_integrate_peak_positions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_peak_positions.md),
[`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md),
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_meta_get()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get.md),
[`nmr_meta_get_column()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get_column.md),
[`nmr_ppm_resolution()`](https://sipss.github.io/AlpsNMR/reference/nmr_ppm_resolution.md),
[`print.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_1D.md)

Other nmr_dataset_peak_table functions:
[`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md),
[`nmr_meta_get()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get.md),
[`nmr_meta_get_column()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get_column.md)

## Examples

``` r
# Load a demo dataset with four samples:
dataset <- system.file("dataset-demo", package = "AlpsNMR")
nmr_dataset <- nmr_read_samples_dir(dataset)

# At first we just have the NMRExperiment column
nmr_meta_get(nmr_dataset, groups = "external")
#> # A tibble: 3 × 1
#>   NMRExperiment
#>   <chr>        
#> 1 10           
#> 2 20           
#> 3 30           
# Get a table with NMRExperiment -> SubjectID
dummy_metadata <- system.file("dataset-demo", "dummy_metadata.xlsx", package = "AlpsNMR")
NMRExp_SubjID <- readxl::read_excel(dummy_metadata, sheet = 1)

NMRExp_SubjID
#> # A tibble: 3 × 3
#>   NMRExperiment SubjectID TimePoint
#>   <chr>         <chr>     <chr>    
#> 1 10            Ana       baseline 
#> 2 20            Ana       3 months 
#> 3 30            Elia      baseline 
# We can link the SubjectID column of the first excel into the dataset
nmr_dataset <- nmr_meta_add(nmr_dataset, NMRExp_SubjID, by = "NMRExperiment")
nmr_meta_get(nmr_dataset, groups = "external")
#> # A tibble: 3 × 3
#>   NMRExperiment SubjectID TimePoint
#>   <chr>         <chr>     <chr>    
#> 1 10            Ana       baseline 
#> 2 20            Ana       3 months 
#> 3 30            Elia      baseline 
# The second excel can use the SubjectID:
SubjID_Age <- readxl::read_excel(dummy_metadata, sheet = 2)
SubjID_Age
#> # A tibble: 2 × 2
#>   SubjectID   Age
#>   <chr>     <dbl>
#> 1 Ana          29
#> 2 Elia          0
# Add the metadata by its SubjectID:
nmr_dataset <- nmr_meta_add(nmr_dataset, SubjID_Age, by = "SubjectID")
# The final loaded metadata:
nmr_meta_get(nmr_dataset, groups = "external")
#> # A tibble: 3 × 4
#>   NMRExperiment SubjectID TimePoint   Age
#>   <chr>         <chr>     <chr>     <dbl>
#> 1 10            Ana       baseline     29
#> 2 20            Ana       3 months     29
#> 3 30            Elia      baseline      0

# Read a tidy excel file:

dataset <- system.file("dataset-demo", package = "AlpsNMR")
nmr_dataset <- nmr_read_samples_dir(dataset)

# At first we just have the NMRExperiment column
nmr_meta_get(nmr_dataset, groups = "external")
#> # A tibble: 3 × 1
#>   NMRExperiment
#>   <chr>        
#> 1 10           
#> 2 20           
#> 3 30           
# Get a table with NMRExperiment -> SubjectID
dummy_metadata <- system.file("dataset-demo", "dummy_metadata.xlsx", package = "AlpsNMR")

nmr_dataset <- nmr_meta_add_tidy_excel(nmr_dataset, dummy_metadata)
# Updated Metadata:
nmr_meta_get(nmr_dataset, groups = "external")
#> # A tibble: 3 × 4
#>   NMRExperiment SubjectID TimePoint   Age
#>   <chr>         <chr>     <chr>     <dbl>
#> 1 10            Ana       baseline     29
#> 2 20            Ana       3 months     29
#> 3 30            Elia      baseline      0
```
