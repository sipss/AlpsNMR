# Plot vip scores of bootstrap

Plot vip scores of bootstrap

## Usage

``` r
plot_vip_scores(vip_means, error, n_samples, plot = TRUE)
```

## Arguments

- vip_means:

  vips means values of bootstraps

- error:

  error tolerated, calculated in the bootstrap

- n_samples:

  number of training samples the bootstrap models were fit on. It sets
  the importance threshold line at `qt(0.975, df = n_samples - 1)`,
  following Afanador, Tran & Buydens (2013), where the cut-off quantile
  `t_{1-alpha/2,n-1}` uses n = number of samples (not the number of
  bootstraps).

- plot:

  A boolean that indicate if results are plotted or not

## Value

A plot of the results or a ggplot object

## Examples

``` r
# Data analysis for a table of integrated peaks

## Generate an artificial nmr_dataset_peak_table:
### Generate artificial metadata:
num_samples <- 64 # use an even number in this example
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

## We will use bootstrap and permutation method for VIPs selection
## in a a k-fold cross validation
# bp_results <- bp_kfold_VIP_analysis(peak_table, # Data to be analized
#                           y_column = "Condition", # Label
#                           k = 3,
#                           ncomp = 1,
#                           nbootstrap = 10)

# message("Selected VIPs are: ", bp_results$importarn_vips)

# plot_vip_scores(bp_results$kfold_results[[1]]$vip_means,
#                bp_results$kfold_results[[1]]$error[1],
#                n_samples = length(bp_results$kfold_index[[1]]))
```
