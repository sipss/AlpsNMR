# Read processed Bruker NMR data

Read processed Bruker NMR data

## Usage

``` r
read_bruker_pdata(
  sample_path,
  pdata_file = NULL,
  pdata_path = "pdata/1",
  all_components = FALSE,
  read_pdata_title = TRUE
)
```

## Arguments

- sample_path:

  A character path of the sample directory

- pdata_file:

  File name of the binary NMR data to load. Usually "1r". If `NULL`, it
  is autodetected based on the dimension

- pdata_path:

  Path from `sample_path` to the preprocessed data

- all_components:

  If `FALSE` load only the real component. Otherwise load the real and
  imaginary components

- read_pdata_title:

  If `TRUE` also reads metadata from pdata title file.

## Value

A list with Bruker NMR processed data
