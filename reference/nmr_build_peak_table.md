# Build a peak table from the clustered peak list

Build a peak table from the clustered peak list

## Usage

``` r
nmr_build_peak_table(peak_data, dataset = NULL)
```

## Arguments

- peak_data:

  A peak list, with the cluster column

- dataset:

  A
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object, to get the metadata

## Value

An
[nmr_dataset_peak_table](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_peak_table.md),
containing the peak table and the annotations

## Examples

``` r
peak_data <- data.frame(
    NMRExperiment = c("10", "10", "20", "20"),
    peak_id = paste0("Peak", 1:4),
    ppm = c(1, 2, 1.1, 2.1),
    area = c(10, 20, 12, 22)
)
clustering_result <- nmr_peak_clustering(peak_data, num_clusters = 2)
peak_data <- clustering_result$peak_data
peak_table <- nmr_build_peak_table(peak_data)
stopifnot(ncol(peak_table) == 2)
```
