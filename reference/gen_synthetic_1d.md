# Generate a synthetic 1D signal with known ground truth

Builds a synthetic spectrum (a smooth baseline plus a number of peaks
plus noise) with a known true baseline and known individual peak areas,
useful for testing baseline-estimation code or trying out
[`psalsa()`](https://sipss.github.io/AlpsNMR/reference/psalsa.md)
without real data. This is what
[`tune_psalsa()`](https://sipss.github.io/AlpsNMR/reference/tune_psalsa.md)
uses internally to build its synthetic tuning pool, matched to the
characteristics of a real example spectrum; called directly, it is a
standalone generator.

## Usage

``` r
gen_synthetic_1d(
  n = 1000,
  density = 0.02,
  fwhm_range = c(10, 30),
  csnr = 0.03,
  peak_shape = c("gaussian", "gex"),
  A = 1,
  seed = 1,
  cap_density = TRUE,
  min_spacing_mult = 2
)
```

## Arguments

- n:

  Number of points.

- density:

  Peaks per point; the number of peaks is `max(3, round(density * n))`
  (further reduced by `cap_density`).

- fwhm_range:

  `c(lo, hi)`: each peak's full width at half maximum (in points) is
  drawn uniformly from this range.

- csnr:

  Noise scale, roughly the noise standard deviation as a fraction of the
  signal amplitude `A`. `0` gives an exactly noise-free signal.

- peak_shape:

  Either `"gaussian"` or `"gex"` (an asymmetric, exponentially-modified
  peak shape typical of chromatography).

- A:

  Amplitude scale for both the baseline and the peaks.

- seed:

  Random seed, for reproducibility.

- cap_density:

  If `TRUE` (the default), reduces the peak count so peaks have room to
  stay non-overlapping (`density` and `fwhm_range` themselves are left
  untouched – only how many peaks are actually placed). Peak centers are
  still placed uniformly at random, so this lowers the odds of overlap
  rather than guaranteeing it; set to `FALSE` to allow arbitrarily
  dense, overlapping peaks.

- min_spacing_mult:

  Only used when `cap_density = TRUE`: the peak count is capped so the
  average spacing between peaks is at least `min_spacing_mult` times the
  widest FWHM in `fwhm_range`.

## Value

A list with:

- `y`:

  the synthetic signal, `baseline + peaks + noise`.

- `baseline`, `peaks`, `noise`:

  the individual components of `y`.

- `sigma`:

  the (possibly position-dependent) noise standard deviation used to
  generate `noise`.

- `peak_info`:

  a data frame with one row per individual peak: `height`, `fwhm`,
  `center`, the index window `lo`:`hi` where that peak alone exceeds
  0.1% of its own height, and `area` (its exact true area,
  [`sum()`](https://rdrr.io/r/base/sum.html) of that peak's own values
  over the full signal).

## Examples

``` r
sig <- AlpsNMR:::gen_synthetic_1d(n = 500, seed = 1)
plot(sig$y, type = "l")
lines(sig$baseline, col = "red")

result <- AlpsNMR:::psalsa(sig$y, lambda = 1e6)
lines(result$baseline, col = "blue", lty = 2)

```
