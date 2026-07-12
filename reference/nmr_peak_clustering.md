# Peak clustering

Peak clustering

## Usage

``` r
nmr_peak_clustering(
  peak_data,
  peak2peak_dist = NULL,
  num_clusters = NULL,
  max_dist_thresh_ppb = NULL,
  verbose = FALSE
)
```

## Arguments

- peak_data:

  The peak list

- peak2peak_dist:

  The distances obtained with
  [nmr_get_peak_distances](https://sipss.github.io/AlpsNMR/reference/nmr_get_peak_distances.md).
  If NULL it is computed from `peak_data`

- num_clusters:

  If you want to fix the number of clusters. Leave `NULL` if you want to
  estimate it

- max_dist_thresh_ppb:

  To estimate the number of clusters, we enforce a limit on how far two
  peaks of the same cluster may be. By default this threshold will be
  computed as 3 times the median peak width (gamma), as given in the
  peak list.

- verbose:

  A logical vector to print additional information

## Value

A list including:

- The `peak_data` with an additional "cluster" column

- cluster: the hierarchical cluster

- num_clusters: an estimation of the number of clusters

- num_cluster_estimation: A list with tables and plots to justify the
  number of cluster estimation

## Examples

``` r
peak_data <- data.frame(
    NMRExperiment = c("10", "10", "20", "20"),
    peak_id = paste0("Peak", 1:4),
    ppm = c(1, 2, 1.1, 2.2),
    gamma_ppb = 100
)
clustering_result <- nmr_peak_clustering(peak_data)
peak_data <- clustering_result$peak_data
stopifnot("cluster" %in% colnames(peak_data))
```
