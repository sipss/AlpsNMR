# Extract parts of an nmr_dataset_peak_table

Extract parts of an nmr_dataset_peak_table

## Usage

``` r
# S3 method for class 'nmr_dataset_peak_table'
x[i]
```

## Arguments

- x:

  an
  [nmr_dataset_peak_table](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_peak_table.md)
  object

- i:

  indices of the samples to keep

## Value

an nmr_dataset_peak_table with the extracted samples

## See also

Other subsetting functions:
[`[.nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset.md),
[`[.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset_1D.md),
[`filter.nmr_dataset_family()`](https://sipss.github.io/AlpsNMR/reference/filter.nmr_dataset_family.md),
[`nmr_pca_outliers_filter()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_filter.md)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
meta <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
metadata <- readxl::read_excel(meta, sheet = 1)
dataset_1D <- nmr_meta_add(dataset_1D, metadata = metadata, by = "NMRExperiment")
metadata <- list(external = dataset_1D[["metadata"]][["external"]])
peak_table <- nmr_data(dataset_1D)
new <- new_nmr_dataset_peak_table(peak_table, metadata)
new[0]
#> An nmr_dataset_peak_table (0 samples, and 45653 peaks) 
```
