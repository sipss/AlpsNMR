# NMR file lister

The function lists samples from the chosen folder required to import and
create a
[nmr_dataset_1D](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
object. The function is based on the
[`fs::dir_ls()`](https://fs.r-lib.org/reference/dir_ls.html) function.

## Usage

``` r
file_lister(dataset_path_nmr, glob)
```

## Arguments

- dataset_path_nmr:

  A character vector of the path where samples are.

- glob:

  A wildcard or globbing pattern common for the samples to be read, for
  example ending with \*0 (spectra acquired by a NOESY sequence often
  end by 0: 10, 20, 30...) or \*s (for example, samples from the
  tutorial in this package) passed on to
  [`grep()`](https://rdrr.io/r/base/grep.html) to filter paths.

## Value

lists of samples from the chosen folder

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
lists_of_samples <- file_lister(dir_to_demo_dataset, "*0")
```
