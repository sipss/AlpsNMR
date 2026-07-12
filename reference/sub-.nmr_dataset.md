# Extract parts of an nmr_dataset

Extract parts of an nmr_dataset

## Usage

``` r
# S3 method for class 'nmr_dataset'
x[i]
```

## Arguments

- x:

  an
  [nmr_dataset](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md)
  object

- i:

  indices of the samples to keep

## Value

an nmr_dataset with the extracted samples

## See also

Other subsetting functions:
[`[.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset_1D.md),
[`[.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset_peak_table.md),
[`filter.nmr_dataset_family()`](https://sipss.github.io/AlpsNMR/reference/filter.nmr_dataset_family.md),
[`nmr_pca_outliers_filter()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_filter.md)

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
dataset2 <- dataset[1:3] # get the first 3 samples
```
