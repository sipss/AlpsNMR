# Bootstrap and permutation over PLS-VIP

Bootstrap and permutation over PLS-VIP on AlpsNMR can be performed on
both
[nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
full spectra as well as
[nmr_dataset_peak_table](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_peak_table.md)
peak tables.

## Usage

``` r
bp_VIP_analysis(
  dataset,
  train_index,
  y_column,
  identity_column = NULL,
  ncomp,
  nbootstrap = 300
)
```

## Arguments

- dataset:

  An
  [nmr_dataset_family](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_family.md)
  object

- train_index:

  set of index used to generate the bootstrap datasets

- y_column:

  A string with the name of the y column (present in the metadata of the
  dataset)

- identity_column:

  `NULL` or a string with the name of the identity column (present in
  the metadata of the dataset). When given, bootstrap resamples are
  drawn by resampling whole `identity_column` groups (e.g. subjects)
  with replacement, rather than individual rows, so that every repeated
  measurement of a resampled subject is kept together, and a multilevel
  (repeated-measures) `plsda` model is fitted (see
  [mixOmics::plsda](https://rdrr.io/pkg/mixOmics/man/plsda.html)'s
  `multilevel` argument), matching what
  [`nmr_data_analysis()`](https://sipss.github.io/AlpsNMR/reference/nmr_data_analysis.md)
  does for the same `identity_column`.

- ncomp:

  number of components used in the plsda models

- nbootstrap:

  number of bootstrap dataset

## Value

A list with the following elements:

- `important_vips`: A list with the important vips selected

- `relevant_vips`: List of vips with some relevance

- `pls_vip`: Pls-VIPs of every bootstrap

- `pls_vip_perm`: Pls-VIPs of every bootstrap with permuted variables

- `pls_vip_means`: Pls-VIPs normaliced differences means

- `pls_vip_score_diff`: Differences of `pls_vip` and `pls_vip_perm`

- `pls_models`: pls models of the diferent bootstraps

- `pls_perm_models`: pls permuted models of the diferent bootstraps

- `classif_rate`: classification rate of the bootstrap models

- `general_model`: pls model trained with all train data

- `general_CR`: classification rate of the `general_model`

- `vips_model`: pls model trained with vips selection over all train
  data

- `vips_CR`: classification rate of the `vips_model`

- `error`: error spected in a t distribution

- `lower_bound`: lower bound of the confidence interval

- `upper_bound`: upper bound of the confidence interval

## Details

Use of the bootstrap and permutation methods for a more robust variable
importance in the projection metric for partial least squares regression

## Examples

``` r
# Data analysis for a table of integrated peaks

## Generate an artificial nmr_dataset_peak_table:
### Generate artificial metadata:
num_samples <- 32 # use an even number in this example
num_peaks <- 20
metadata <- data.frame(
    NMRExperiment = as.character(1:num_samples),
    Condition = rep(c("A", "B"), times = num_samples / 2)
)

### The matrix with peaks
peak_means <- runif(n = num_peaks, min = 300, max = 600)
peak_sd <- runif(n = num_peaks, min = 30, max = 60)
peak_matrix <- mapply(function(mu, sd) rnorm(num_samples, mu, sd),
    mu = peak_means, sd = peak_sd
)
colnames(peak_matrix) <- paste0("Peak", 1:num_peaks)

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

## We will use a double cross validation, splitting the samples with random
## subsampling both in the external and internal validation.
## The classification model will be a PLSDA, exploring at maximum 3 latent
## variables.
## The best model will be selected based on the area under the ROC curve
methodology <- plsda_auroc_vip_method(ncomp = 3)
model <- nmr_data_analysis(
    peak_table,
    y_column = "Condition",
    identity_column = NULL,
    external_val = list(iterations = 1, test_size = 0.25),
    internal_val = list(iterations = 3, test_size = 0.25),
    data_analysis_method = methodology
)
## Area under ROC for each outer cross-validation iteration:
model$outer_cv_results_digested$auroc
#> # A tibble: 1 × 3
#>   cv_outer_iteration ncomp   auc
#>                <int> <int> <dbl>
#> 1                  1     1     1

## The number of components for the bootstrap models is selected
ncomps <- model$outer_cv_results$`1`$model$ncomp
train_index <- model$train_test_partitions$outer$`1`$outer_train

# Bootstrap and permutation for VIP selection
bp_VIPS <- bp_VIP_analysis(peak_table, # Data to be analyzed
    train_index,
    y_column = "Condition",
    ncomp = ncomps,
    nbootstrap = 10
)
#> Warning: No VIPs are ranked as important
#> ℹ Try increasing the number of bootstrap iterations
#> Warning: bp_VIP_analysis: no relevant_vips found
#> ℹ You may try increasing the number of bootstraps
```
