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
#' @param peak_shape One of `"lorentzian"` (the default; the peak shape
#'   typical of NMR spectra, the same Cauchy/Lorentzian used by
#'   [peaklist_fit_lorentzians()]), `"gaussian"`, or `"gex"` (an asymmetric,
#'   exponentially-modified peak shape typical of chromatography).
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
#'       peak alone exceeds 0.1% of its own height, `area` (its exact
#'       true area, `sum()` of that peak's own values over the full signal),
#'       and `a`/`b` (the `"gex"` shape parameters, `NA` for other shapes) --
#'       enough to exactly reconstruct that one peak's own curve later,
#'       without storing a dense per-peak matrix.}
#'     \item{`peak_shape`}{the `peak_shape` used, needed alongside
#'       `peak_info` to reconstruct individual peaks.}
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
                              peak_shape = c("lorentzian", "gaussian", "gex"),
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
    a_j <- NA_real_; b_j <- NA_real_
    pk_j <- if (peak_shape == "gaussian") {
      gauss_peak_1d(x, centers[j], fw, h[j])
    } else if (peak_shape == "lorentzian") {
      lorentz_peak_1d(x, centers[j], fw, h[j])
    } else {
      a_j <- stats::runif(1, 0.5, 2); b_j <- stats::runif(1, 5, 8)
      gex_peak_1d(x, centers[j] - fw / 2, centers[j] + fw / 2, h[j], a_j, b_j)
    }
    peaks <- peaks + pk_j
    above <- which(pk_j > 1e-3 * h[j])
    lo <- if (length(above)) min(above) else max(1, round(centers[j]))
    hi <- if (length(above)) max(above) else min(n, round(centers[j]))
    peak_rows[[j]] <- data.frame(idx = j, height = h[j], fwhm = fw, center = centers[j],
                                  lo = lo, hi = hi, area = sum(pk_j), a = a_j, b = b_j)
  }
  peak_info <- do.call(rbind, peak_rows)

  sigma <- A * csnr * (0.6 + 0.4 * xf)
  noise <- if (csnr == 0) numeric(n) else stats::rnorm(n, 0, 1) * sigma

  list(y = baseline + peaks + noise, baseline = baseline, peaks = peaks, noise = noise,
       sigma = sigma, peak_info = peak_info, peak_shape = peak_shape)
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

#' Generate a synthetic 1D signal whose peak density/width vary by region
#'
#' Like [gen_synthetic_1d()], but instead of one constant `density` and
#' `fwhm_range` for the whole signal, takes a `region_profile` (as produced by
#' `pool_signal_stats_regions_1d()` in `psalsa_tune.R`) so different stretches
#' of the synthetic spectrum can be crowded or sparse, matching how a real
#' spectrum's peak density actually varies from region to region (e.g. a
#' crowded aliphatic region vs. a sparse downfield region) -- something a
#' single pooled density/fwhm blends away. This is what [tune_psalsa()] uses
#' internally when called with `num_regions` set; [gen_synthetic_1d()] (no
#' region awareness) remains the default.
#'
#' @param n Number of points.
#' @param region_profile A data frame with one row per region.: `frac_lo`,
#'   `frac_hi` (the region's span as a fraction of the signal, `0..1`),
#'   `density` (peaks per point within the region), `fwhm_q1`, `fwhm_q3` (that
#'   region's FWHM range, in points). A region with `density <= 0` or `NA`/
#'   non-positive `fwhm_q1`/`fwhm_q3` (no peaks ever observed there) gets no
#'   synthetic peaks placed in it -- an empty stretch is real information, not
#'   missing data. An optional `amplitude` column gives that region's own
#'   typical peak height scale (used in place of `A` for that region's peaks
#'   only -- see `A` below); if absent, or `NA`/non-positive for a given
#'   region, that region's peaks fall back to `A` -- byte-identical to not
#'   having this column at all, so existing callers are unaffected.
#' @param csnr,peak_shape,A,seed,cap_density,min_spacing_mult As in
#'   [gen_synthetic_1d()], applied per-region (`min_spacing_mult` against that
#'   region's own `fwhm_q3`, not the global widest FWHM). `A` and `csnr` are
#'   used AS GIVEN (not overridden by `region_profile$amplitude`) for the
#'   smooth baseline component and the noise level (`sigma <- A * csnr *
#'   ...`) -- keeping those tied to the signal-wide `A`/`csnr` rather than a
#'   region's own amplitude is what keeps the ABSOLUTE noise level constant
#'   across regions of differing peak height, matching real
#'   (instrument/electronic) noise, which doesn't scale with local peak
#'   amplitude the way `A`-scaled peak heights should.
#'
#' @return As in [gen_synthetic_1d()].
#' @noRd
gen_synthetic_1d_regions <- function(n, region_profile, csnr = 0.03,
                                      peak_shape = c("lorentzian", "gaussian", "gex"),
                                      A = 1, seed = 1, cap_density = TRUE, min_spacing_mult = 2) {
  peak_shape <- match.arg(peak_shape)
  set.seed(seed)
  x <- seq_len(n); xf <- x / n
  baseline <- A * (1 + 0.5 * xf + 0.3 * xf^2 + 0.15 * sin(2 * pi * 2.3 * xf + stats::runif(1, -pi, pi)))

  peaks <- numeric(n); peak_rows <- list()
  for (r in seq_len(nrow(region_profile))) {
    reg <- region_profile[r, ]
    lo <- max(1L, floor(reg$frac_lo * n) + 1L)
    hi <- min(n, floor(reg$frac_hi * n))
    reg_n <- hi - lo + 1L
    if (reg_n < 1 || is.na(reg$density) || reg$density <= 0 ||
      is.na(reg$fwhm_q1) || is.na(reg$fwhm_q3) || reg$fwhm_q3 <= 0) {
      next
    }

    npk_r <- round(reg$density * reg_n)
    if (cap_density) {
      max_npk_r <- floor(reg_n / (min_spacing_mult * reg$fwhm_q3))
      npk_r <- min(npk_r, max_npk_r)
    }
    if (npk_r < 1) next

    reg_amplitude <- if (!is.null(region_profile$amplitude) &&
      !is.na(region_profile$amplitude[r]) && region_profile$amplitude[r] > 0) {
      region_profile$amplitude[r]
    } else {
      A
    }
    h_r <- stats::rlnorm(npk_r, meanlog = log(0.35 * reg_amplitude), sdlog = log(6))
    centers_r <- stats::runif(npk_r, lo, hi)
    for (j in seq_len(npk_r)) {
      fw <- stats::runif(1, reg$fwhm_q1, reg$fwhm_q3)
      a_j <- NA_real_; b_j <- NA_real_
      pk_j <- if (peak_shape == "gaussian") {
        gauss_peak_1d(x, centers_r[j], fw, h_r[j])
      } else if (peak_shape == "lorentzian") {
        lorentz_peak_1d(x, centers_r[j], fw, h_r[j])
      } else {
        a_j <- stats::runif(1, 0.5, 2); b_j <- stats::runif(1, 5, 8)
        gex_peak_1d(x, centers_r[j] - fw / 2, centers_r[j] + fw / 2, h_r[j], a_j, b_j)
      }
      peaks <- peaks + pk_j
      above <- which(pk_j > 1e-3 * h_r[j])
      lo_j <- if (length(above)) min(above) else max(1, round(centers_r[j]))
      hi_j <- if (length(above)) max(above) else min(n, round(centers_r[j]))
      peak_rows[[length(peak_rows) + 1]] <- data.frame(
        idx = length(peak_rows) + 1, height = h_r[j], fwhm = fw, center = centers_r[j],
        lo = lo_j, hi = hi_j, area = sum(pk_j), a = a_j, b = b_j
      )
    }
  }
  peak_info <- if (length(peak_rows)) {
    do.call(rbind, peak_rows)
  } else {
    data.frame(
      idx = integer(0), height = numeric(0), fwhm = numeric(0), center = numeric(0),
      lo = integer(0), hi = integer(0), area = numeric(0), a = numeric(0), b = numeric(0)
    )
  }

  sigma <- A * csnr * (0.6 + 0.4 * xf)
  noise <- if (csnr == 0) numeric(n) else stats::rnorm(n, 0, 1) * sigma

  list(
    y = baseline + peaks + noise, baseline = baseline, peaks = peaks, noise = noise,
    sigma = sigma, peak_info = peak_info, peak_shape = peak_shape
  )
}

#' Reconstruct one synthetic peak's own curve from its `peak_info` row
#'
#' Exactly reproduces one peak's own contribution (no baseline, no noise, no
#' other peaks) over `x`, using the same shape function [gen_synthetic_1d()]
#' used to build it in the first place. Used to fractionally attribute a
#' baseline-corrected *overlapping* region between the peaks that share it,
#' rather than requiring peaks to be isolated to score their area recovery
#' at all -- see `peak_area_errors_1d()` in `psalsa_tune.R`.
#'
#' @param x Points to evaluate the peak at (typically just its own `lo:hi`
#'   window, not the full signal).
#' @param peak_shape As in [gen_synthetic_1d()].
#' @param height,fwhm,center As in one row of `peak_info`.
#' @param a,b The `"gex"` shape parameters (ignored for other shapes).
#' @return A numeric vector the same length as `x`.
#' @noRd
reconstruct_peak_1d <- function(x, peak_shape, height, fwhm, center, a = NA_real_, b = NA_real_) {
  if (peak_shape == "gaussian") {
    gauss_peak_1d(x, center, fwhm, height)
  } else if (peak_shape == "lorentzian") {
    lorentz_peak_1d(x, center, fwhm, height)
  } else {
    gex_peak_1d(x, center - fwhm / 2, center + fwhm / 2, height, a, b)
  }
}
