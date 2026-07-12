# K-fold bootstrap and permutation over PLS-VIP

Bootstrap and permutation over PLS-VIP on AlpsNMR can be performed on
both
[nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
full spectra as well as
[nmr_dataset_peak_table](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_peak_table.md)
peak tables.

## Usage

``` r
bp_kfold_VIP_analysis(dataset, y_column, k = 4, ncomp = 3, nbootstrap = 300)
```

## Arguments

- dataset:

  An
  [nmr_dataset_family](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_family.md)
  object

- y_column:

  A string with the name of the y column (present in the metadata of the
  dataset)

- k:

  Number of folds, recomended between 4 to 10

- ncomp:

  number of components for the bootstrap models

- nbootstrap:

  number of bootstrap dataset

## Value

A list with the following elements:

- `important_vips`: A list with the important vips selected

- `relevant_vips`: List of vips with some relevance

- `wilcoxon_vips`: List of vips that pass a wilcoxon test

- `vip_means`: Means of the vips scores

- `vip_score_plot`: plot of the vips scores

- `kfold_resuls`: results of the k
  [bp_VIP_analysis](https://sipss.github.io/AlpsNMR/reference/bp_VIP_analysis.md)

- `kfold_index`: list of index of partitions of the folds

## Details

Use of the bootstrap and permutation methods for a more robust variable
importance in the projection metric for partial least squares
regression, in a k-fold cross validation

## Examples

``` r
# Data analysis for a table of integrated peaks
set.seed(42)
## Generate an artificial nmr_dataset_peak_table:
### Generate artificial metadata:
num_samples <- 64 # use an even number in this example
num_peaks <- 10
metadata <- data.frame(
    NMRExperiment = as.character(1:num_samples),
    Condition = sample(rep(c("A", "B"), times = num_samples / 2), num_samples)
)

### The matrix with peaks
peak_means <- runif(n = num_peaks, min = 300, max = 600)
peak_sd <- runif(n = num_peaks, min = 30, max = 60)
peak_matrix <- mapply(function(mu, sd) rnorm(num_samples, mu, sd),
    mu = peak_means, sd = peak_sd
)
colnames(peak_matrix) <- paste0("Peak", 1:num_peaks)
rownames(peak_matrix) <- paste0("Sample", 1:num_samples)

## Artificial differences depending on the condition:
peak_matrix[metadata$Condition == "A", "Peak2"] <-
    peak_matrix[metadata$Condition == "A", "Peak2"] + 70

peak_matrix[metadata$Condition == "A", "Peak6"] <-
    peak_matrix[metadata$Condition == "A", "Peak6"] - 60

### The nmr_dataset_peak_table
peak_table <- new_nmr_dataset_peak_table(
    peak_table = peak_matrix,
    metadata = list(external = metadata)
)

## We will use bootstrap and permutation method for VIPs selection
## in a a k-fold cross validation
bp_results <- bp_kfold_VIP_analysis(peak_table, # Data to be analyzed
    y_column = "Condition", # Label
    k = 2,
    ncomp = 1,
    nbootstrap = 5
)
#> Warning: No VIPs are ranked as important
#> ℹ Try increasing the number of bootstrap iterations

message("Selected VIPs are: ", bp_results$important_vips)
#> Selected VIPs are: 
```
