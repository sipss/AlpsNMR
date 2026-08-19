# Tune PSALSA parameters from one or more example spectra

[`psalsa()`](https://sipss.github.io/AlpsNMR/reference/psalsa.md) needs
`lambda`, `p` and `k` to be chosen for the data at hand, and good values
depend on the signal's noise level, peak width and peak density.
`tune_psalsa()` picks them automatically from example spectra, without
requiring a known baseline: it analyzes the noise level and a rough peak
width/density summary of `y`, generates a batch of *synthetic* spectra
with matching characteristics (for which the true baseline and peak
areas are known by construction, since they were generated), and
searches for the `lambda`/`p`/`k` that best recovers those synthetic
peaks' areas. The resulting parameters are then applied to `y` to
produce the returned baseline.

## Usage

``` r
tune_psalsa(
  y,
  peak_shape = c("gaussian", "gex"),
  n_synthetic = 10,
  optim_maxit = 150,
  optim_reltol = 1e-06
)
```

## Arguments

- y:

  A numeric vector with one example spectrum, or a list of numeric
  vectors with several example spectra (their characteristics are pooled
  together before tuning). Lists of differing lengths are supported.

- peak_shape:

  Either `"gaussian"` or `"gex"` (an asymmetric, exponentially-modified
  peak shape typical of chromatography). Peak shape is often known from
  the instrument/technique used to acquire `y`, so it is left as an
  explicit argument rather than inferred from the data.

- n_synthetic:

  Number of synthetic spectra generated for the search.

- optim_maxit:

  Maximum number of
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) (Nelder-Mead)
  iterations.

- optim_reltol:

  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) relative
  convergence tolerance.

## Value

If `y` is a single numeric vector, a list with:

- `baseline`, `corrected`:

  as in
  [`psalsa()`](https://sipss.github.io/AlpsNMR/reference/psalsa.md), for
  `y` with the tuned parameters.

- `lambda`, `p`, `k`:

  the tuned parameters, reusable via
  `psalsa(other_y, lambda = lambda, p = p, k = k)` on similar spectra
  without tuning again.

If `y` is a list, `baseline` and `corrected` are lists (one element per
input spectrum, using the same tuned `lambda`/`p`/`k` for all of them).

## Details

The search balances recovery across small, medium and large peaks (by
true area) so that a few large peaks don't dominate the objective at the
expense of small ones; it falls back to a plain baseline-accuracy
objective when the synthetic batch doesn't contain enough well-separated
peaks to define those three groups (e.g. very dense or very wide peaks,
which tend to overlap).

## See also

[`psalsa()`](https://sipss.github.io/AlpsNMR/reference/psalsa.md), the
function being tuned.

## Examples

``` r
x <- seq_len(300)
baseline <- 10 + 5 * sin(x / 60)
peaks <- 60 * exp(-((x - 90)^2) / (2 * 4^2)) +
  25 * exp(-((x - 220)^2) / (2 * 6^2))
set.seed(1)
y <- baseline + peaks + rnorm(length(x), 0, 0.5)

result <- AlpsNMR:::tune_psalsa(y)
plot(y, type = "l")
lines(result$baseline, col = "red")

```
