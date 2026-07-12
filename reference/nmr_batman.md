# Batman helpers

Batman helpers

## Usage

``` r
nmr_batman_write_options(
  bopts,
  batman_dir = "BatmanInput",
  filename = "batmanOptions.txt"
)

nmr_batman_export_dataset(
  nmr_dataset,
  batman_dir = "BatmanInput",
  filename = "NMRdata.txt"
)

nmr_batman_multi_data_user_hmdb(
  batman_dir = "BatmanInput",
  filename = "multi_data_user.csv"
)

nmr_batman_multi_data_user(
  multiplet_table,
  batman_dir = "BatmanInput",
  filename = "multi_data_user.csv"
)

nmr_batman_metabolites_list(
  metabolite_names,
  batman_dir = "BatmanInput",
  filename = "metabolitesList.csv"
)
```

## Arguments

- bopts:

  Batman options

- batman_dir:

  Batman input directorye

- filename:

  Filename to use, inside `batman_dir`

- nmr_dataset:

  An
  [nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  object

- multiplet_table:

  A data frame, like the
  [hmdb](https://sipss.github.io/AlpsNMR/reference/hmdb.md) dataset

- metabolite_names:

  A character vector of the metabolite names to consider

## Value

These are helper functions to make Batman tests easier

## See also

Other batman functions:
[`nmr_batman_options()`](https://sipss.github.io/AlpsNMR/reference/nmr_batman_options.md)

## Examples

``` r
bopts <- nmr_batman_options()
# nmr_batman_write_options(bopts)

dataset_1D <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
# nmr_batman_export_dataset(dataset_1D)

message("Use of multi_data_user_hmdb")
#> Use of multi_data_user_hmdb
# multi_data_user_hmdb <- nmr_batman_multi_data_user_hmdb()
hmdb <- NULL
# utils::data("hmdb", package = "AlpsNMR", envir = environment())
# hmdb <- nmr_batman_multi_data_user(hmbd)

metabolite_names <- c("alanine", "glucose")
# metabolite_names <- nmr_batman_metabolites_list(metabolite_names)
```
