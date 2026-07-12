# Fit lorentzians to each peak to estimate areas

The different methods are available for benchmarking while developing,
we should pick one.

## Usage

``` r
peaklist_fit_lorentzians(
  peak_data,
  nmr_dataset,
  amplitude_method = c("intensity", "2nd_derivative", "intensity_without_baseline"),
  refine_peak_model = c("none", "peak", "2nd_derivative")
)
```

## Arguments

- peak_data:

  The peak data

- nmr_dataset:

  The nmr_dataset object with the data. This function for now assumes
  nmr_dataset is NOT be baseline corrected

- amplitude_method:

  The method to estimate the amplitude. It may be:

  - `"intensity"`. The amplitude of the peak is proportional to the raw
    intensity at the apex. This is a bad estimation if the intensity
    includes a baseline, because the amplitude of the peak will be
    overestimated

  - `"2nd_derivative"`: The amplitude of the peak is proportional to the
    second derivative of the raw intensity signal at the apex. This
    method aims to correct the "intensity" method, since it is expected
    that the baseline will be mostly removed when considering the 2nd
    derivative of the spectrum. The 2nd derivative is calculated with a
    2nd order Savitzky-Golay filter of 21 points.

  - `"intensity_without_baseline"`: A baseline is estimated on the whole
    spectra and subtracted from it. Then the peak amplitude is
    proportional to the corrected intensity at the apex (as in the
    "intensity" method).

- refine_peak_model:

  Whether a non linear least squares fitting should be used to refine
  the estimated parameters. It can be:

  - `"none"`: Do not refine using nls.

  - `"peak"`: Use a lorentzian peak model and the baseline corrected
    spectra.

  - `"2nd_derivative"`:

## Value

The given data frame `peak_data`, with added columns:

- inflection points,

- gamma

- area

- a norm_rmse fitting error

As well as some attributes

- "errors": A data frame with any error in the peak fitting

- "fit_baseline": Whether the method used has any consideration for the
  baseline of the signal (maybe not very useful attribute)

- "method_description": A textual description of what we did, to include
  it in plots

## Details

- gamma is estimated using the inflection points of the signal and
  fitting them to the lorentzian inflection points

- \$A\$ is estimated using the `amplitude_method` below

- The peak position (\$x_0\$) is given in `peak_data`

Those estimations may be refined with non-linear least squares using
`refine_peak_model`. If the nls does not converge, the initial
estimations are kept. Convergence -and other nls errors- are saved for
further reference and diagnostic. Use `attr(peak_data_fitted, "errors")`
to retreive the error messages, where `peak_data_fitted` is assumed to
be the output of this function. The refining improves gamma, \$A\$ and
\$x_0\$.

The baseline estimation (when calculated, see the arguments) is set to
Asymmetric Least Squares with lambda = 6, p=0.05, maxit=20 and it is
probably not optimal... yet.
