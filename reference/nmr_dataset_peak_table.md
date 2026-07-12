# nmr_dataset_peak_table (S3 class)

An `nmr_dataset_peak_table` represents a peak table with metadata. It is
defined as an S3 class, and it can be treated as a regular list.

## Usage

``` r
# S3 method for class 'nmr_dataset_peak_table'
as.data.frame(x, ...)
```

## Arguments

- x:

  An nmr_dataset_peak_table object,

- ...:

  ignored

## Value

A data frame with the sample metadata and the peak table

## Details

- `metadata`: A list of data frames. Each data frame contains metadata.
  Usually the list only has one data frame named "external".

- `peak_table`: A matrix with one sample on each row and the peaks in
  the columns

## Methods (by generic)

- `as.data.frame(nmr_dataset_peak_table)`: Convert to a data frame
