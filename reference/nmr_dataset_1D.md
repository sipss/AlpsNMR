# nmr_dataset_1D (S3 class)

An `nmr_dataset_1D` represents a set of 1D interpolated NMR samples. It
is defined as an S3 class, and it can be treated as a regular list.

## Details

It currently has the following elements:

- `metadata`: A list of data frames. Each data frame contains metadata
  of a given area (acquisition parameters, preprocessing parameters,
  general sample information...)

- `axis`: A numeric vector with the chemical shift axis in ppm.

- `data_1r`: A matrix with one sample on each row and the chemical
  shifts in the columns.

## Examples

``` r
# Create a random spectra matrix
nsamp <- 12
npoints <- 20
dummy_ppm_axis <- seq(from = 0.2, to = 10, length.out = npoints)
dummy_spectra_matrix <- matrix(runif(nsamp * npoints), nrow = nsamp, ncol = npoints)
metadata <- list(external = data.frame(
    NMRExperiment = paste0("Sample", 1:12),
    DummyClass = c("a", "b")
))
dummy_nmr_dataset_1D <- new_nmr_dataset_1D(
    ppm_axis = dummy_ppm_axis,
    data_1r = dummy_spectra_matrix,
    metadata = metadata
)
```
