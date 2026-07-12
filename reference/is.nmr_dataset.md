# Object is of [nmr_dataset](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md) class

Object is of
[nmr_dataset](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md)
class

## Usage

``` r
is.nmr_dataset(x)
```

## Arguments

- x:

  An object

## Value

`TRUE` if the object is an
[nmr_dataset](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md),
`FALSE` otherwise

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
is(dataset)
#> [1] "nmr_dataset"
```
