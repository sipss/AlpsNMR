# Peak list: Create an `accepted` column based on some criteria

Peak list: Create an `accepted` column based on some criteria

## Usage

``` r
peaklist_accept_peaks(
  peak_data,
  nmr_dataset,
  nrmse_max = Inf,
  area_min = 0,
  area_max = Inf,
  ppm_min = -Inf,
  ppm_max = Inf,
  keep_rejected = TRUE,
  verbose = FALSE
)
```

## Arguments

- peak_data:

  The peak list (a data frame)

- nmr_dataset:

  The nmr_dataset where the peak_data was computed from

- nrmse_max:

  The normalized root mean squared error of the lorentzian peak fitting
  must be less than or equal to this value

- area_min:

  Peak areas must be larger or equal to this value

- area_max:

  Peak areas must be smaller or equal to this value

- ppm_min:

  The peak apex must be above this value

- ppm_max:

  The peak apex must be below this value

- keep_rejected:

  If `FALSE`, removes those peaks that do not satisfy the criteria and
  remove the accepted column (since all would be accepted)

- verbose:

  Print informational message

## Value

The `peak_data`, with a new `accepted` column (or maybe some filtered
rows)

## Examples

``` r
# Fake data:
nmr_dataset <- new_nmr_dataset_1D(
    1:10,
    matrix(c(1:5, 4:2, 3, 0), nrow = 1),
    list(external = data.frame(NMRExperiment = "10"))
)
peak_data <- data.frame(
    peak_id = c("Peak1", "Peak2"),
    NMRExperiment = c("10", "10"),
    ppm = c(5, 9),
    pos = c(5, 9),
    intensity = c(5, 3),
    ppm_infl_min = c(3, 8),
    ppm_infl_max = c(7, 10),
    gamma_ppb = c(1, 1),
    area = c(25, 3),
    norm_rmse = c(0.01, 0.8)
)
# Create the accepted column:
peak_data <- peaklist_accept_peaks(peak_data, nmr_dataset, area_min = 10, keep_rejected = FALSE)
stopifnot(identical(peak_data$peak_id, "Peak1"))
```
