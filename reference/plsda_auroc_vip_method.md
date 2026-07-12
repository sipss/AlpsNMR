# Method for nmr_data_analysis (PLSDA model with AUROC and VIP outputs)

Method for nmr_data_analysis (PLSDA model with AUROC and VIP outputs)

## Usage

``` r
plsda_auroc_vip_method(ncomp, auc_increment_threshold = 0.05)
```

## Arguments

- ncomp:

  Max. number of latent variables to explore in the PLSDA analysis

- auc_increment_threshold:

  Choose the number of latent variables when the AUC does not increment
  more than this threshold.

## Value

Returns an object to be used with
[nmr_data_analysis](https://sipss.github.io/AlpsNMR/reference/nmr_data_analysis.md)
to perform a (optionally multilevel) PLS-DA model, using the area under
the ROC curve as figure of merit to determine the optimum number of
latent variables.

## Examples

``` r
method <- plsda_auroc_vip_method(3)
```
