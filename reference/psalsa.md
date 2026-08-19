# Baseline estimation with Peaked Signal's Asymmetric Least Squares

Estimates the baseline of one or more spectra using the PSALSA (Peaked
Signal's Asymmetric Least Squares Algorithm) method. PSALSA is a variant
of asymmetric least squares that limits the influence of intense peaks
on the estimated baseline, which improves the result for signals
containing sharp, high peaks.

## Usage

``` r
psalsa(
  spectra,
  lambda = 1e+07,
  p = 0.001,
  k = -1,
  maxit = 25,
  k_epsilon = 1e-06
)
```

## Arguments

- spectra:

  Either a numeric vector containing one spectrum, or a matrix
  containing one spectrum per row.

- lambda:

  Smoothing parameter (generally `1e5` - `1e8`).

- p:

  Asymmetry parameter.

- k:

  Peak height parameter, controlling how strongly the weights of points
  above the baseline decay with their height. Usually about 5% of the
  maximum intensity. When `k = -1` (the default) it is set to one
  twentieth of the maximum intensity of each spectrum.

- maxit:

  Maximum number of iterations.

- k_epsilon:

  Floors `k` at `k_epsilon * diff(range(spectra))`, so an estimated or
  supplied `k` that is (near) zero can't make the asymmetric weight
  `exp(-d/k)` collapse to 0 for essentially any positive residual –
  which happens for a genuinely near-noise-free signal, since `k` is
  usually derived from an estimated noise level (see
  [`tune_psalsa()`](https://sipss.github.io/AlpsNMR/reference/tune_psalsa.md)),
  silently degenerating the baseline toward the signal's minimum. Scaled
  by the data's own range (rather than a fixed absolute value) so the
  floor stays proportionate regardless of the signal's amplitude.

## Value

A list with two elements, each with the same dimensions as `spectra`:

- `baseline`:

  the estimated baseline.

- `corrected`:

  the baseline-corrected signal, i.e. `spectra - baseline`.

## Details

Like ordinary asymmetric least squares, PSALSA assigns different weights
to the points above and below an iteratively estimated baseline: the
asymmetry parameter `p` (`0 <= p <= 1`) is the weight for points below
the baseline, whereas points above it normally receive weight `1 - p`.
The difference is that for points above the baseline the weight decays
exponentially with the height of the point above the current estimate,
controlled by `k`. This prevents intense peaks from pulling the baseline
upwards. The parameter `lambda` controls the amount of smoothing: the
larger it is, the smoother the baseline will be. Iteration stops once
the set of points above the baseline no longer changes, or after `maxit`
iterations.

The baseline is fitted with a weighted second-order Whittaker smoother.

## References

Oller-Moreno, S., Pardo, A., Jimenez-Soto, J. M., Samitier, J., Marco,
S. (2014). "Adaptive Asymmetric Least Squares baseline estimation for
analytical instruments". 2014 IEEE 11th International Multi-Conference
on Systems, Signals & Devices (SSD14), 1-5.
[doi:10.1109/SSD.2014.6808837](https://doi.org/10.1109/SSD.2014.6808837)

## Examples

``` r
x <- seq_len(200)
baseline <- 10 + 5 * sin(x / 40)
peak <- 80 * exp(-((x - 120)^2) / (2 * 4^2))
y <- baseline + peak

result <- AlpsNMR:::psalsa(y)
plot(y, type = "l")
lines(result$baseline, col = "red")

```
