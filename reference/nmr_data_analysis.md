# Data analysis

Data analysis on AlpsNMR can be performed on both
[nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
full spectra as well as
[nmr_dataset_peak_table](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_peak_table.md)
peak tables.

## Usage

``` r
nmr_data_analysis(
  dataset,
  y_column,
  identity_column,
  external_val,
  internal_val,
  data_analysis_method,
  .enable_parallel = TRUE
)
```

## Arguments

- dataset:

  An
  [nmr_dataset_family](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_family.md)
  object

- y_column:

  A string with the name of the y column (present in the metadata of the
  dataset)

- identity_column:

  `NULL` or a string with the name of the identity column (present in
  the metadata of the dataset).

- external_val, internal_val:

  A list with two elements: `iterations` and `test_size`. See
  [random_subsampling](https://sipss.github.io/AlpsNMR/reference/random_subsampling.md)
  for further details

- data_analysis_method:

  An
  [nmr_data_analysis_method](https://sipss.github.io/AlpsNMR/reference/nmr_data_analysis_method.md)
  object

- .enable_parallel:

  Set to `FALSE` to disable parallellization.

## Value

A list with the following elements:

- `train_test_partitions`: A list with the indices used in train and
  test on each of the cross-validation iterations

- `inner_cv_results`: The output returned by `train_evaluate_model` on
  each inner cross-validation

- `inner_cv_results_digested`: The output returned by
  `choose_best_inner`.

- `outer_cv_results`: The output returned by `train_evaluate_model` on
  each outer cross-validation

- `outer_cv_results_digested`: The output returned by
  `train_evaluate_model_digest_outer`.

## Details

The workflow consists of a double cross validation strategy using random
subsampling for splitting into train and test sets. The classification
model and the metric to choose the best model can be customized (see
[`new_nmr_data_analysis_method()`](https://sipss.github.io/AlpsNMR/reference/nmr_data_analysis_method.md)),
but for now only a PLSDA classification model with a best area under ROC
curve metric is implemented (see the examples here and
[plsda_auroc_vip_method](https://sipss.github.io/AlpsNMR/reference/plsda_auroc_vip_method.md))

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
    external_val = list(iterations = 3, test_size = 0.25),
    internal_val = list(iterations = 3, test_size = 0.25),
    data_analysis_method = methodology
)
## Area under ROC for each outer cross-validation iteration:
model$outer_cv_results_digested$auroc
#> # A tibble: 3 × 3
#>   cv_outer_iteration ncomp   auc
#>                <int> <int> <dbl>
#> 1                  1     1 0.733
#> 2                  2     1 0.5  
#> 3                  3     1 0.6  
## Rank Product of the Variable Importance in the Projection
## (Lower means more important)
sort(model$outer_cv_results_digested$vip_rankproducts)
#>     Peak2     Peak6    Peak16     Peak1    Peak11     Peak7    Peak13    Peak18 
#>  2.000000  2.080084  2.620741  3.107233  4.308869  5.768998  7.958114  9.491220 
#>    Peak14     Peak8    Peak12     Peak3    Peak20     Peak9    Peak17    Peak19 
#>  9.813199 11.292432 11.617571 12.324480 12.557072 13.782348 13.959064 14.851875 
#>     Peak4     Peak5    Peak10    Peak15 
#> 15.036946 15.740609 16.108636 17.621736 
```
