#' Generate a synthetic 1D signal with known ground truth
#'
#' Builds a synthetic spectrum (a smooth baseline plus a number of peaks plus
#' noise) with a known true baseline and known individual peak areas, useful
#' for testing baseline-estimation code or trying out [psalsa()] without real
#' data. This is what [tune_psalsa()] uses internally to build its synthetic
#' tuning pool, matched to the characteristics of a real example spectrum;
#' called directly, it is a standalone generator.
#'
#' @param n Number of points.
#' @param density Peaks per point; the number of peaks is
#'   `max(3, round(density * n))` (further reduced by `cap_density`).
#' @param fwhm_range `c(lo, hi)`: each peak's full width at half maximum (in
#'   points) is drawn uniformly from this range.
#' @param csnr Noise scale, roughly the noise standard deviation as a
#'   fraction of the signal amplitude `A`. `0` gives an exactly noise-free
#'   signal.
#' @param peak_shape One of `"gaussian"`, `"gex"` (an asymmetric,
#'   exponentially-modified peak shape typical of chromatography), or
#'   `"lorentzian"` (the peak shape typical of NMR spectra; the same
#'   Cauchy/Lorentzian used by [peaklist_fit_lorentzians()]).
#' @param A Amplitude scale for both the baseline and the peaks.
#' @param seed Random seed, for reproducibility.
#' @param cap_density If `TRUE` (the default), reduces the peak count so
#'   peaks have room to stay non-overlapping (`density` and `fwhm_range`
#'   themselves are left untouched -- only how many peaks are actually
#'   placed). Peak centers are still placed uniformly at random, so this
#'   lowers the odds of overlap rather than guaranteeing it; set to `FALSE`
#'   to allow arbitrarily dense, overlapping peaks.
#' @param min_spacing_mult Only used when `cap_density = TRUE`: the peak
#'   count is capped so the average spacing between peaks is at least
#'   `min_spacing_mult` times the widest FWHM in `fwhm_range`.
#'
#' @return A list with:
#'   \describe{
#'     \item{`y`}{the synthetic signal, `baseline + peaks + noise`.}
#'     \item{`baseline`, `peaks`, `noise`}{the individual components of `y`.}
#'     \item{`sigma`}{the (possibly position-dependent) noise standard
#'       deviation used to generate `noise`.}
#'     \item{`peak_info`}{a data frame with one row per individual peak:
#'       `height`, `fwhm`, `center`, the index window `lo`:`hi` where that
#'       peak alone exceeds 0.1% of its own height, and `area` (its exact
#'       true area, `sum()` of that peak's own values over the full signal).}
#'   }
#'
#'
#' @examples
#' sig <- AlpsNMR:::gen_synthetic_1d(n = 500, seed = 1)
#' plot(sig$y, type = "l")
#' lines(sig$baseline, col = "red")
#'
#' result <- AlpsNMR:::psalsa(sig$y, lambda = 1e6)
#' lines(result$baseline, col = "blue", lty = 2)
#'
gen_synthetic_1d <- function(n = 1000, density = 0.02, fwhm_range = c(10, 30), csnr = 0.03,
                              peak_shape = c("gaussian", "gex", "lorentzian"),
                              A = 1, seed = 1, cap_density = TRUE, min_spacing_mult = 2) {
  peak_shape <- match.arg(peak_shape)
  set.seed(seed)
  x <- seq_len(n); xf <- x / n
  baseline <- A * (1 + 0.5 * xf + 0.3 * xf^2 + 0.15 * sin(2 * pi * 2.3 * xf + stats::runif(1, -pi, pi)))

  npk <- max(3, round(density * n))
  if (cap_density) {
    max_npk <- max(3, floor(n / (min_spacing_mult * fwhm_range[2])))
    npk <- min(npk, max_npk)
  }
  h <- stats::rlnorm(npk, meanlog = log(0.35 * A), sdlog = log(6))
  centers <- stats::runif(npk, 0.02 * n, 0.98 * n)

  peaks <- numeric(n); peak_rows <- vector("list", npk)
  for (j in seq_len(npk)) {
    fw <- stats::runif(1, fwhm_range[1], fwhm_range[2])
    pk_j <- if (peak_shape == "gaussian") {
      gauss_peak_1d(x, centers[j], fw, h[j])
    } else if (peak_shape == "lorentzian") {
      lorentz_peak_1d(x, centers[j], fw, h[j])
    } else {
      a <- stats::runif(1, 0.5, 2); b <- stats::runif(1, 5, 8)
      gex_peak_1d(x, centers[j] - fw / 2, centers[j] + fw / 2, h[j], a, b)
    }
    peaks <- peaks + pk_j
    above <- which(pk_j > 1e-3 * h[j])
    lo <- if (length(above)) min(above) else max(1, round(centers[j]))
    hi <- if (length(above)) max(above) else min(n, round(centers[j]))
    peak_rows[[j]] <- data.frame(idx = j, height = h[j], fwhm = fw, center = centers[j],
                                  lo = lo, hi = hi, area = sum(pk_j))
  }
  peak_info <- do.call(rbind, peak_rows)

  sigma <- A * csnr * (0.6 + 0.4 * xf)
  noise <- if (csnr == 0) numeric(n) else stats::rnorm(n, 0, 1) * sigma

  list(y = baseline + peaks + noise, baseline = baseline, peaks = peaks, noise = noise,
       sigma = sigma, peak_info = peak_info)
}

## Exponentially-modified peak shape (asymmetric rise, exponential tail),
## typical of chromatographic peaks.
#' @noRd
gex_peak_1d <- function(x, t0, tm, h, a, b) {
  u <- (x - t0) / (tm - t0)
  pk <- numeric(length(x)); pos <- u > 0
  pk[pos] <- h * u[pos]^(b - 1) * exp((b - 1) / a * (1 - u[pos]^a))
  pk[!is.finite(pk)] <- 0
  pk
}

#' @noRd
gauss_peak_1d <- function(x, center, fwhm, h) {
  sigma <- fwhm / (2 * sqrt(2 * log(2)))
  h * exp(-((x - center)^2) / (2 * sigma^2))
}

## Height/fwhm-parametrized wrapper around the package's canonical
## lorentzian() (area_estimation.R, also used by peaklist_fit_lorentzians()):
## gamma there is the half-width-at-half-maximum, so fwhm = 2*gamma, and A is
## chosen so the peak's value at its center equals h (A/(pi*gamma) = h).
#' @noRd
lorentz_peak_1d <- function(x, center, fwhm, h) {
  gamma <- fwhm / 2
  lorentzian(x, x0 = center, gamma = gamma, A = h * pi * gamma)
}
