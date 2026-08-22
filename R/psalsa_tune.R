#' Tune PSALSA parameters from one or more example spectra
#'
#' `psalsa()` needs `lambda`, `p` and `k` to be chosen for the data at hand, and
#' good values depend on the signal's noise level, peak width and peak density.
#' `tune_psalsa()` picks them automatically from example spectra, without
#' requiring a known baseline: it analyzes the noise level and a rough peak
#' width/density summary of `y`, generates a batch of *synthetic* spectra with
#' matching characteristics (for which the true baseline and peak areas are
#' known by construction, since they were generated), and searches for the
#' `lambda`/`p`/`k` that best recovers those synthetic peaks' areas. The
#' resulting parameters are then applied to `y` to produce the returned
#' baseline.
#'
#' The search balances recovery across small, medium and large peaks (by true
#' area) so that a few large peaks don't dominate the objective at the expense
#' of small ones; it falls back to a plain baseline-accuracy objective when
#' the synthetic batch doesn't contain enough well-separated peaks to define
#' those three groups (e.g. very dense or very wide peaks, which tend to
#' overlap).
#'
#' @param y A numeric vector with one example spectrum, or a list of numeric
#'   vectors with several example spectra (their characteristics are pooled
#'   together before tuning). Lists of differing lengths are supported.
#' @param peak_shape One of `"lorentzian"` (the default; the peak shape
#'   typical of NMR spectra), `"gaussian"`, or `"gex"` (an asymmetric,
#'   exponentially-modified peak shape typical of chromatography). Peak shape
#'   is often known from the instrument/technique used to acquire `y`, so it
#'   is left as an explicit argument rather than inferred from the data.
#' @param n_synthetic Number of synthetic spectra generated for the search.
#' @param optim_maxit Maximum number of [stats::optim()] (Nelder-Mead)
#'   iterations.
#' @param optim_reltol [stats::optim()] relative convergence tolerance.
#' @param num_regions If set, characterizes peak density/width separately in
#'   this many contiguous regions of `y` (instead of one signal-wide, pooled
#'   density/width) and generates the synthetic tuning pool to match that
#'   region-by-region variation -- useful when `y`'s peak density is
#'   noticeably uneven across the spectrum (e.g. crowded aliphatic vs. sparse
#'   downfield regions), since a single blended density otherwise makes every
#'   synthetic region look equally busy. `lambda`/`p`/`k` are still tuned as
#'   single, signal-wide values either way; this only changes how realistic
#'   the synthetic pool is. `NULL` (the default) keeps the original
#'   signal-wide characterization.
#'
#' @return If `y` is a single numeric vector, a list with:
#'   \describe{
#'     \item{`baseline`, `corrected`}{as in [psalsa()], for `y` with the tuned
#'       parameters.}
#'     \item{`lambda`, `p`, `k`}{the tuned parameters, reusable via
#'       `psalsa(other_y, lambda = lambda, p = p, k = k)` on similar spectra
#'       without tuning again.}
#'   }
#'   If `y` is a list, `baseline` and `corrected` are lists (one element per
#'   input spectrum, using the same tuned `lambda`/`p`/`k` for all of them).
#'
#' @seealso [psalsa()], the function being tuned.
#'
#' @examples
#' x <- seq_len(300)
#' baseline <- 10 + 5 * sin(x / 60)
#' peaks <- 60 * exp(-((x - 90)^2) / (2 * 4^2)) +
#'   25 * exp(-((x - 220)^2) / (2 * 6^2))
#' set.seed(1)
#' y <- baseline + peaks + rnorm(length(x), 0, 0.5)
#'
#' result <- AlpsNMR:::tune_psalsa(y)
#' plot(y, type = "l")
#' lines(result$baseline, col = "red")
#'
tune_psalsa <- function(y, peak_shape = c("lorentzian", "gaussian", "gex"), n_synthetic = 10,
                         optim_maxit = 150, optim_reltol = 1e-6, num_regions = NULL) {
  peak_shape <- match.arg(peak_shape)
  y_list <- if (is.list(y)) y else list(y)
  if (is.null(num_regions)) {
    pooled_stats <- pool_signal_stats_1d(y_list)
    pool <- generate_synthetic_pool_1d(pooled_stats, y_list, n_synthetic = n_synthetic, peak_shape = peak_shape)
  } else {
    pooled_stats <- pool_signal_stats_regions_1d(y_list, num_regions = num_regions)
    pool <- generate_synthetic_pool_1d_regions(pooled_stats, y_list, n_synthetic = n_synthetic, peak_shape = peak_shape)
  }
  tuned <- tune_psalsa_params_1d(pool, optim_maxit = optim_maxit, optim_reltol = optim_reltol)
  k <- tuned$k_mult * pooled_stats$noise_sd

  if (is.list(y)) {
    fits <- lapply(y, function(yi) psalsa(yi, lambda = tuned$lambda, p = tuned$p, k = k))
    baseline <- lapply(fits, `[[`, "baseline")
    corrected <- lapply(fits, `[[`, "corrected")
  } else {
    fit <- psalsa(y, lambda = tuned$lambda, p = tuned$p, k = k)
    baseline <- fit$baseline
    corrected <- fit$corrected
  }
  list(baseline = baseline, corrected = corrected, lambda = tuned$lambda, p = tuned$p, k = k)
}

#' Tune position-varying PSALSA parameters from one or more example spectra
#'
#' Like [tune_psalsa()], but instead of a single, signal-wide `lambda`/`p`/`k`,
#' tunes each of `num_regions` regions independently (via
#' `tune_psalsa_region_params_1d()`: each region gets its own synthetic pool,
#' matched to that region's own peak density/width, and its own
#' `tune_psalsa_params_1d()` search) and combines the per-region optima into
#' smooth, position-varying `lambda`/`p`/`k` profiles (one value per point),
#' via a natural cubic spline over each parameter's own transformed scale --
#' see `spatial_profile_1d()`. [psalsa()] is then run ONCE per spectrum with
#' these profiles, using the position-varying `diff2_penalty_weighted()`
#' generalization of the smoothing penalty.
#'
#' Useful when a single, signal-wide `lambda`/`p`/`k` (as tuned by
#' [tune_psalsa()], even with `num_regions` set there merely to build a more
#' realistic synthetic pool) still isn't enough -- e.g. because one region's
#' own baseline behaviour genuinely needs a different smoothing/peak-
#' resistance strength than another, not just a more realistic peak density
#' during tuning.
#'
#' Each region's own search is ridge-regularized toward a shared anchor (via
#' `tune_psalsa_params_1d()`'s own regularization mechanism, at a much
#' stronger weight -- see `region_lambda_k_weight`/`region_p_weight`)
#' computed from a single signal-wide tune of the WHOLE spectrum first --
#' "the pack" -- rather than toward `tune_psalsa_params_1d()`'s generic
#' literature-default prior. A region with few peaks of its own has a flat,
#' weakly-constrained objective; verified empirically that anchoring alone,
#' at the SAME regularization weight the whole-spectrum search uses, is not
#' enough to stop such a region's knot from still collapsing to the same
#' "disable peak protection" degenerate optimum described in
#' `tune_psalsa_params_1d()`'s own comments (`lambda` collapsing toward 0,
#' `k` exploding) -- a region has far less data to resist that escape hatch
#' than the whole spectrum does. The stronger per-region weights keep that
#' region close to the pack unless its own data genuinely justifies moving
#' away from it.
#'
#' A region with no observed peaks contributes no knot (there's nothing local
#' to tune against); the spline fills that stretch in smoothly from its
#' neighbouring regions' knots instead. If fewer than 2 regions have enough
#' signal to tune at all, this falls back to plain, signal-wide [tune_psalsa()].
#'
#' @inheritParams tune_psalsa
#' @param num_regions Number of contiguous regions to characterize and tune
#'   independently.
#' @param p_max Upper bound used for `p`'s logit-space spline interpolation
#'   (see `spatial_profile_1d()`); should match `tune_psalsa_params_1d()`'s own
#'   `p_max` (its default, `0.05`, is used here too).
#' @param region_lambda_k_weight,region_p_weight Ridge-regularization weights
#'   for each region's own search toward the signal-wide anchor -- the same
#'   role as `tune_psalsa_params_1d()`'s own `lambda_k_weight`/`p_weight`, but
#'   defaulting much stronger (verified necessary; see Details) since a single
#'   region typically has far less data to constrain the search than the
#'   whole spectrum does.
#' @param min_peaks Before tuning, adjacent regions with fewer than this many
#'   estimated peaks are merged together (see `merge_sparse_regions_1d()`),
#'   so every region actually tuned has enough real peaks of its own to
#'   constrain the search -- addressing sparse-region collapse at its root
#'   cause (too little data) rather than relying only on
#'   `region_lambda_k_weight`/`region_p_weight` to paper over it.
#'
#' @return If `y` is a single numeric vector, a list with:
#'   \describe{
#'     \item{`baseline`, `corrected`}{as in [psalsa()], for `y` with the tuned
#'       position-varying parameters.}
#'     \item{`lambda`, `p`, `k`}{the tuned position-varying parameter profiles,
#'       each the same length as `y`, reusable via `psalsa(other_y, lambda =
#'       lambda, p = p, k = k)` on similarly-lengthed spectra without tuning
#'       again (resample with `resample_profile_1d()` first if the length
#'       differs).}
#'     \item{`region_params`}{a data frame with the per-region tuned values
#'       the profiles were spline-interpolated from (`region`, `frac_mid`,
#'       `lambda`, `p`, `k_mult`), useful for inspecting what the tuning
#'       actually found region by region.}
#'     \item{`noise_sd`}{the signal-wide noise level `k_mult` was scaled by
#'       to get `k`. `region_params` + `noise_sd` are enough to rebuild
#'       `lambda`/`p`/`k` later without re-tuning -- much smaller to store
#'       than the profiles themselves -- via `psalsa_region_params_to_profile()`.}
#'   }
#'   If `y` is a list, `baseline` and `corrected` are lists (one element per
#'   input spectrum, using the same tuned profiles -- resampled to that
#'   spectrum's own length where needed).
#'
#' @seealso [tune_psalsa()], [psalsa()].
#'
#' @examples
#' x <- seq_len(600)
#' baseline <- 10 + 5 * sin(x / 100)
#' # A crowded region (many small peaks) and a sparse region (one big peak):
#' peaks <- Reduce(`+`, lapply(seq(50, 250, by = 20), function(c) {
#'   8 * exp(-((x - c)^2) / (2 * 3^2))
#' })) + 60 * exp(-((x - 500)^2) / (2 * 6^2))
#' set.seed(1)
#' y <- baseline + peaks + rnorm(length(x), 0, 0.5)
#'
#' result <- AlpsNMR:::tune_psalsa_spatial(y, num_regions = 6, n_synthetic = 5)
#' plot(y, type = "l")
#' lines(result$baseline, col = "red")
#'
tune_psalsa_spatial <- function(y, peak_shape = c("lorentzian", "gaussian", "gex"), n_synthetic = 10,
                                 optim_maxit = 150, optim_reltol = 1e-6, num_regions = 20, p_max = 0.05,
                                 region_lambda_k_weight = 2, region_p_weight = 10, min_peaks = 15) {
  peak_shape <- match.arg(peak_shape)
  y_list <- if (is.list(y)) y else list(y)
  n <- round(stats::median(vapply(y_list, length, integer(1))))

  ## Signal-wide anchor ("the pack"): a single, ordinary tune_psalsa() pass
  ## over the WHOLE spectrum, used only to regularize each region's own
  ## search below (see tune_psalsa_region_params_1d()'s comments) -- not
  ## itself returned or used as a baseline fit.
  pooled_global <- pool_signal_stats_1d(y_list)
  pool_global <- generate_synthetic_pool_1d(pooled_global, y_list, n_synthetic = n_synthetic, peak_shape = peak_shape)
  tuned_global <- tune_psalsa_params_1d(pool_global, p_max = p_max, optim_maxit = optim_maxit, optim_reltol = optim_reltol)
  theta0_anchor <- c(log(tuned_global$lambda), stats::qlogis(tuned_global$p / p_max), log(tuned_global$k_mult))

  pooled_regions <- pool_signal_stats_regions_1d(y_list, num_regions = num_regions)
  region_params <- tune_psalsa_region_params_1d(
    pooled_regions, n_total = n, peak_shape = peak_shape,
    n_synthetic = n_synthetic, optim_maxit = optim_maxit, optim_reltol = optim_reltol,
    theta0 = theta0_anchor, p_max = p_max,
    region_lambda_k_weight = region_lambda_k_weight, region_p_weight = region_p_weight,
    min_peaks = min_peaks
  )

  if (is.null(region_params) || nrow(region_params) < 2) {
    ## Too few regions had enough signal to tune independently (e.g. a very
    ## sparse or very short spectrum) -- a spline needs at least 2 knots, so
    ## fall back to a single, signal-wide tune rather than a degenerate one.
    return(tune_psalsa(
      y, peak_shape = peak_shape, n_synthetic = n_synthetic,
      optim_maxit = optim_maxit, optim_reltol = optim_reltol
    ))
  }

  profiles <- psalsa_region_params_to_profile(n, region_params, noise_sd = pooled_regions$noise_sd, p_max = p_max)
  lambda_profile <- profiles$lambda
  p_profile <- profiles$p
  k_profile <- profiles$k

  fit_one <- function(yi) {
    ni <- length(yi)
    psalsa(
      yi,
      lambda = resample_profile_1d(lambda_profile, ni),
      p = resample_profile_1d(p_profile, ni),
      k = resample_profile_1d(k_profile, ni)
    )
  }

  if (is.list(y)) {
    fits <- lapply(y, fit_one)
    baseline <- lapply(fits, `[[`, "baseline")
    corrected <- lapply(fits, `[[`, "corrected")
  } else {
    fit <- fit_one(y)
    baseline <- fit$baseline
    corrected <- fit$corrected
  }

  list(
    baseline = baseline, corrected = corrected,
    lambda = lambda_profile, p = p_profile, k = k_profile,
    region_params = region_params, noise_sd = pooled_regions$noise_sd
  )
}

#' Rebuild position-varying PSALSA profiles from a `tune_psalsa_spatial()` region table
#'
#' `tune_psalsa_spatial()`'s tuned `lambda`/`p`/`k` profiles are as long as
#' the spectrum they were tuned on (impractical to store verbatim, e.g. as a
#' literal in a script or vignette for later reuse without re-tuning), but
#' they are entirely determined by its much smaller `region_params` table
#' (one row per region -- typically a few dozen at most) together with
#' `noise_sd`. This rebuilds the full-length profiles from just those two
#' small, easily-stored pieces, via the same spline interpolation
#' `tune_psalsa_spatial()` itself uses (see `spatial_profile_1d()`) --
#' `psalsa_region_params_to_profile(n, tuned$region_params, tuned$noise_sd)`
#' reproduces `list(lambda = tuned$lambda, p = tuned$p, k = tuned$k)` exactly
#' for `n` equal to the length `tuned` was originally computed for.
#'
#' @param n Length of the profile to build (the length of the spectra it will
#'   be applied to via `psalsa()`).
#' @param region_params A data frame with one row per region: `frac_mid`
#'   (that region's midpoint, as a fraction 0-1 of the spectrum),
#'   `lambda`, `p`, `k_mult` -- exactly the `region_params` element of a
#'   `tune_psalsa_spatial()` result.
#' @param noise_sd The `noise_sd` element of the same `tune_psalsa_spatial()`
#'   result (used to turn `k_mult` back into an absolute `k`).
#' @param p_max Upper bound used for `p`'s logit-space interpolation; must
#'   match the `p_max` the original `tune_psalsa_spatial()` call used (its
#'   default, `0.05`, is used here too).
#'
#' @return A list with `lambda`, `p`, `k` -- numeric vectors of length `n`,
#'   directly reusable via `psalsa(y, lambda = lambda, p = p, k = k)` (or
#'   `nmr_baseline_estimation()`'s `lambda`/`p`/`k` arguments) without
#'   calling `tune_psalsa_spatial()` again.
#'
#' @examples
#' # A crowded region (many small peaks) and a sparse region (one big peak),
#' # each with enough peaks of its own to tune independently:
#' n <- 800
#' x <- seq_len(n)
#' baseline <- 5 + 2 * sin(x / 200)
#' peak_centers <- c(seq(10, 190, by = 10), seq(210, 390, by = 10))
#' peaks <- Reduce(`+`, lapply(peak_centers, function(ctr) {
#'   15 * exp(-((x - ctr)^2) / (2 * 1.5^2))
#' }))
#' set.seed(1)
#' y <- baseline + peaks + rnorm(n, 0, 0.2)
#'
#' tuned <- AlpsNMR:::tune_psalsa_spatial(y, num_regions = 4, n_synthetic = 5)
#' # tuned$region_params (a small table) and tuned$noise_sd are cheap to
#' # store (e.g. as literals in a script), unlike tuned$lambda/p/k directly:
#' profiles <- AlpsNMR:::psalsa_region_params_to_profile(length(y), tuned$region_params, tuned$noise_sd)
#' stopifnot(all.equal(profiles$lambda, tuned$lambda))
#'
psalsa_region_params_to_profile <- function(n, region_params, noise_sd, p_max = 0.05) {
  lambda <- spatial_profile_1d(n, region_params$frac_mid, region_params$lambda, transform = "log")
  p <- spatial_profile_1d(n, region_params$frac_mid, region_params$p, transform = "logit", p_max = p_max)
  k_mult <- spatial_profile_1d(n, region_params$frac_mid, region_params$k_mult, transform = "log")
  list(lambda = lambda, p = p, k = k_mult * noise_sd)
}

## Robust noise sd from the MAD of 2nd differences. For pure white noise e,
## Var(e[i-1] - 2*e[i] + e[i+1]) = 6*sigma^2, so sigma_hat = MAD(2nd-diff) /
## sqrt(6) (MAD is already scaled to be a consistent sigma estimator under
## normality). A smooth baseline contributes ~0 to 2nd differences and the
## median-based MAD suppresses sparse peaks, so this is robust to both.
#' @noRd
noise_est_1d <- function(y) {
  n <- length(y)
  d2 <- y[1:(n - 2)] - 2 * y[2:(n - 1)] + y[3:n]
  stats::mad(d2) / sqrt(6)
}

## Rolling low-percentile filter: for a series of overlapping windows, take
## the q-th percentile of y within the window as a control point at the
## window's center, then linearly interpolate. Unlike a smoothing pass (which
## only reflects the current baseline estimate's own dynamics), this is
## structurally biased toward sitting below peaks from the very first pass,
## giving informative residuals immediately -- used only as a cheap reference
## for peak detection below, not as a baseline estimate in its own right.
##
## window_size is capped at max_window (not just window_frac*n unbounded):
## a signal's baseline curvature scale is usually set by the
## instrument/technique, NOT by how many samples happen to have been
## recorded, so sizing the window as a pure FRACTION of n makes it grow
## without bound for a longer signal even when the true local structure
## doesn't change scale at all. Verified this breaks peak detection badly:
## on a signal with a fixed ~200-sample structural period, tiling it out to
## 4x length grew the (uncapped) window from 10 to 40 samples -- an
## increasingly poor match to the true local curvature -- and inflated
## detected peaks from 2 (true 1) to 127 (true 4) on the same underlying
## per-period content, corrupting the density/fwhm estimates used to build
## the synthetic tuning pool.
#' @noRd
percentile_init_1d <- function(y, window_frac = 1 / 20, shift_frac = 0.5, q = 0.10, max_window = 20) {
  n <- length(y)
  window_size <- min(max(3, round(n * window_frac)), max_window)
  window_shift <- max(1, round(window_size * shift_frac))
  starts <- seq(1, n, by = window_shift)

  centers <- numeric(0); vals <- numeric(0)
  for (s in starts) {
    e <- min(s + window_size - 1, n)
    if (e < s) next
    idx <- s:e
    centers <- c(centers, mean(idx))
    vals <- c(vals, as.numeric(stats::quantile(y[idx], q, names = FALSE)))
  }
  if (centers[1] > 1) { centers <- c(1, centers); vals <- c(vals[1], vals) }
  if (utils::tail(centers, 1) < n) { centers <- c(centers, n); vals <- c(vals, utils::tail(vals, 1)) }

  stats::approx(centers, vals, xout = seq_len(n), rule = 2)$y
}

## Simple strict-local-max peak detection on the residual against
## percentile_init_1d() -- deliberately simple (no non-max suppression for
## bumpy peak tops), since this only needs a BALLPARK peak count/width, not
## precise peak-picking.
#' @noRd
detect_peaks_1d <- function(y, height_mult = 5) {
  z0 <- percentile_init_1d(y)
  noise_sd <- noise_est_1d(y)
  r <- y - z0
  n <- length(r)
  thresh <- height_mult * noise_sd

  is_max <- logical(n)
  is_max[2:(n - 1)] <- r[2:(n - 1)] > r[1:(n - 2)] & r[2:(n - 1)] > r[3:n]
  candidates <- which(is_max & r > thresh)

  rows <- lapply(candidates, function(i) {
    half <- r[i] / 2
    lo <- i; while (lo > 1 && r[lo] > half) lo <- lo - 1
    hi <- i; while (hi < n && r[hi] > half) hi <- hi + 1
    data.frame(idx = i, height = r[i], fwhm = hi - lo)
  })
  if (length(rows) == 0) return(data.frame(idx = integer(0), height = numeric(0), fwhm = numeric(0)))
  do.call(rbind, rows)
}

#' @noRd
summarize_peaks_1d <- function(peaks, n) {
  if (nrow(peaks) == 0) {
    return(list(n_peaks = 0, density = 0, fwhm_q1 = NA_real_, fwhm_q2 = NA_real_, fwhm_q3 = NA_real_))
  }
  q <- stats::quantile(peaks$fwhm, c(0.25, 0.5, 0.75), names = FALSE)
  list(n_peaks = nrow(peaks), density = nrow(peaks) / n,
       fwhm_q1 = q[1], fwhm_q2 = q[2], fwhm_q3 = q[3])
}

## One real signal -> noise_sd/csnr (continuous signal-to-noise ratio) + a
## rough peak summary, used to generate a matched synthetic pool. A (the
## amplitude scale used to convert noise_sd into a relative csnr) prefers
## abs(median(percentile_init_1d(y))), falling back to a peak-height-based
## estimate and finally to noise_sd itself when the baseline is genuinely
## close to zero (a plain baseline estimate can dip slightly negative on pure
## noise there, which would otherwise give a nonsensical negative csnr).
#' @noRd
analyze_signal_1d <- function(y, height_mult = 5) {
  noise_sd <- noise_est_1d(y)
  z0 <- percentile_init_1d(y)
  peaks <- detect_peaks_1d(y, height_mult = height_mult)
  pk_summary <- summarize_peaks_1d(peaks, length(y))

  A_baseline <- abs(stats::median(z0))
  A_peaks <- if (nrow(peaks) > 0) stats::median(peaks$height) / 0.35 else NA_real_
  A <- if (A_baseline >= noise_sd) A_baseline
       else if (!is.na(A_peaks) && A_peaks >= noise_sd) A_peaks
       else noise_sd
  csnr <- noise_sd / A

  c(list(noise_sd = noise_sd, csnr = csnr), pk_summary)
}

## Pools analyze_signal_1d() across one or more real input signals: stacks
## each signal's stats into one row, takes the column-wise median (robust to
## one outlier signal in the batch). na.rm=TRUE since fwhm_q1/q2/q3 are NA for
## any input with zero detected peaks.
#' @noRd
pool_signal_stats_1d <- function(y_list) {
  stats_list <- lapply(y_list, analyze_signal_1d)
  stats_df <- do.call(rbind, lapply(stats_list, as.data.frame))
  as.list(sapply(stats_df, stats::median, na.rm = TRUE))
}

## TRUE if a peak's [lo,hi] window doesn't overlap any OTHER peak's window in
## the same signal -- direct pairwise check (peak counts are small, at most a
## few dozen per signal, so O(n^2) is trivial and correct; a sort-and-check-
## adjacent-only shortcut is wrong when one peak's window nests around another
## with an unrelated peak's bound sorted in between).
#' @noRd
peak_is_isolated_1d <- function(peak_info) {
  n <- nrow(peak_info)
  if (n < 2) return(rep(TRUE, n))
  iso <- rep(TRUE, n)
  for (i in seq_len(n)) {
    overlaps <- peak_info$lo[i] <= peak_info$hi & peak_info$lo <= peak_info$hi[i]
    overlaps[i] <- FALSE
    if (any(overlaps)) iso[i] <- FALSE
  }
  iso
}

## Peak-area recovery via fractional attribution. A hard lo:hi window
## (summing the baseline-corrected signal over just that peak's own window)
## double-counts/misattributes area wherever peaks overlap, and has to
## exclude overlapping peaks entirely to stay valid -- see
## peak_is_isolated_1d(). Since these are SYNTHETIC peaks, we know each
## peak's own true contribution at every point (peaks_total is the exact sum
## of every peak alone, no baseline/noise), so instead of a hard window we
## split the corrected signal at each point in proportion to how much each
## peak TRULY contributes there:
##   w_i(x) = pk_i(x) / peaks_total(x)
##   est_area_i = sum_x w_i(x) * corrected(x)
## Where a peak is genuinely isolated, w_i is 1 inside its window and 0
## outside, so this reduces exactly to the old sum(corrected[lo:hi]) there --
## strictly more general, not a different method for the isolated case.
## peaks_total(x) >= pk_i(x) > 0 for every x inside that peak's own lo:hi
## (lo:hi is defined as exactly where pk_i itself exceeds 1e-3*height, and
## peaks stack additively with non-negative shapes), so no zero-division
## guard is needed.
#' @noRd
peak_area_errors_1d <- function(y, z_est, peak_info, peaks_total, peak_shape) {
  if (is.null(peak_info) || nrow(peak_info) == 0) return(peak_info)
  corrected <- as.vector(y) - as.vector(z_est)
  peaks_total <- as.vector(peaks_total)
  est_area <- vapply(seq_len(nrow(peak_info)), function(i) {
    idx <- peak_info$lo[i]:peak_info$hi[i]
    pk_i <- reconstruct_peak_1d(
      x = idx, peak_shape = peak_shape,
      height = peak_info$height[i], fwhm = peak_info$fwhm[i], center = peak_info$center[i],
      a = peak_info$a[i], b = peak_info$b[i]
    )
    w_i <- pk_i / peaks_total[idx]
    sum(w_i * corrected[idx])
  }, numeric(1))
  pct_err <- 100 * (est_area - peak_info$area) / peak_info$area
  cbind(peak_info, est_area = est_area, pct_error = pct_err, abs_pct_error = abs(pct_err),
        isolated = peak_is_isolated_1d(peak_info))
}

## Small/medium/large tiers by true peak area, using cut points computed once
## from the pool's WHOLE peak population (not just isolated peaks -- now that
## peak_area_errors_1d() can score overlapping peaks too, there's no reason
## to leave them out of the tier boundaries), so tiers are comparable across
## the different synthetic signals scored by the same tuning run.
#' @noRd
compute_peak_size_breaks_1d <- function(signals) {
  areas <- unlist(lapply(signals, function(s) {
    pi <- s$peak_info
    if (is.null(pi) || nrow(pi) == 0) return(NULL)
    pi$area
  }))
  stats::quantile(areas, c(1 / 3, 2 / 3), na.rm = TRUE)
}

#' @noRd
classify_peak_size_1d <- function(area, breaks) {
  cut(area, breaks = c(-Inf, breaks, Inf), labels = c("small", "medium", "large"))
}

## RMSE (as a percentage of scale) restricted to points where the peak
## contribution is small relative to noise (peaks < 0.5*sigma) -- the
## baseline is directly observable there, so this is a useful fallback
## objective when the pool has too few peaks (of any kind) for size-tiered
## peak-area scoring, or when the size tiers themselves are degenerate.
#' @noRd
floor_rmse_1d <- function(est, truth, peaks, sigma, scale) {
  mask <- as.vector(peaks) < 0.5 * as.vector(sigma)
  if (sum(mask) < 5) mask <- rep(TRUE, length(truth))
  100 * sqrt(mean((as.vector(est)[mask] - as.vector(truth)[mask])^2)) / scale
}

## pooled_stats (pool_signal_stats_1d() output) + the real input signals (for
## their lengths) -> a list of n_synthetic gen_synthetic_1d() draws, one per
## seed 1..n_synthetic, all matched to the pooled scenario.
#' @noRd
generate_synthetic_pool_1d <- function(pooled_stats, y_list, n_synthetic = 10, peak_shape = "gaussian") {
  n <- round(stats::median(vapply(y_list, length, integer(1))))
  lapply(seq_len(n_synthetic), function(seed) {
    gen_synthetic_1d(n = n, density = pooled_stats$density,
                      fwhm_range = c(pooled_stats$fwhm_q1, pooled_stats$fwhm_q3),
                      csnr = pooled_stats$csnr, peak_shape = peak_shape, seed = seed)
  })
}

## n, num_regions -> a data frame with one row per region (region, lo, hi):
## a contiguous, near-equal-sized integer partition of 1:n. num_regions is
## capped at n so a degenerately short signal never gets empty regions.
#' @noRd
compute_region_bounds_1d <- function(n, num_regions) {
  num_regions <- max(1, min(num_regions, n))
  edges <- floor(seq(0, n, length.out = num_regions + 1))
  data.frame(region = seq_len(num_regions), lo = edges[-length(edges)] + 1L, hi = edges[-1])
}

## Like analyze_signal_1d(), but characterizes peak density/fwhm PER REGION
## instead of pooling them into one signal-wide density/fwhm -- a spectrum's
## peak density is rarely uniform (e.g. a crowded aliphatic region vs. a
## sparse downfield region), and a single blended density/fwhm can't
## reproduce that when generating synthetic tuning data. noise_sd/csnr stay
## signal-wide (not per-region): regions from num_regions=20 are typically
## too short for a robust per-region noise_est_1d(), and noise level is
## primarily an instrument/acquisition characteristic, not a per-region one.
## region bounds are stored as FRACTIONS of n (frac_lo/frac_hi), not raw
## indices, so per-region stats stay comparable/poolable across real input
## signals of different lengths.
#' @noRd
analyze_signal_regions_1d <- function(y, num_regions = 20, height_mult = 5) {
  n <- length(y)
  noise_sd <- noise_est_1d(y)
  z0 <- percentile_init_1d(y)
  peaks <- detect_peaks_1d(y, height_mult = height_mult)

  A_baseline <- abs(stats::median(z0))
  A_peaks <- if (nrow(peaks) > 0) stats::median(peaks$height) / 0.35 else NA_real_
  A <- if (A_baseline >= noise_sd) A_baseline
       else if (!is.na(A_peaks) && A_peaks >= noise_sd) A_peaks
       else noise_sd
  csnr <- noise_sd / A

  bounds <- compute_region_bounds_1d(n, num_regions)
  region_rows <- lapply(seq_len(nrow(bounds)), function(r) {
    lo <- bounds$lo[r]; hi <- bounds$hi[r]
    in_region <- peaks$idx >= lo & peaks$idx <= hi
    pk_summary <- summarize_peaks_1d(peaks[in_region, , drop = FALSE], hi - lo + 1)
    data.frame(
      region = r, frac_lo = (lo - 1) / n, frac_hi = hi / n,
      n_peaks = pk_summary$n_peaks, density = pk_summary$density,
      fwhm_q1 = pk_summary$fwhm_q1, fwhm_q2 = pk_summary$fwhm_q2, fwhm_q3 = pk_summary$fwhm_q3
    )
  })
  list(noise_sd = noise_sd, csnr = csnr, regions = do.call(rbind, region_rows))
}

## Pools analyze_signal_regions_1d() across one or more real input signals:
## noise_sd/csnr pool the same way pool_signal_stats_1d() does (column-wise
## median), and each region's density/fwhm pools (median, na.rm=TRUE) against
## the SAME region index across signals -- valid because region bounds are
## fractional, so "region 7 of 20" means the same relative stretch of the
## spectrum regardless of a given signal's own length. A region with zero
## detected peaks in every pooled signal keeps density 0 and fwhm NA (that's
## real information -- see gen_synthetic_1d_regions(), which places no peaks
## there rather than guessing a width for peaks that were never observed).
#' @noRd
pool_signal_stats_regions_1d <- function(y_list, num_regions = 20, height_mult = 5) {
  per_signal <- lapply(y_list, analyze_signal_regions_1d, num_regions = num_regions, height_mult = height_mult)
  noise_sd <- stats::median(vapply(per_signal, `[[`, numeric(1), "noise_sd"))
  csnr <- stats::median(vapply(per_signal, `[[`, numeric(1), "csnr"))

  regions_stacked <- do.call(rbind, lapply(per_signal, `[[`, "regions"))
  pooled_regions <- do.call(rbind, lapply(split(regions_stacked, regions_stacked$region), function(sub) {
    data.frame(
      region = sub$region[1], frac_lo = stats::median(sub$frac_lo), frac_hi = stats::median(sub$frac_hi),
      density = stats::median(sub$density, na.rm = TRUE),
      fwhm_q1 = stats::median(sub$fwhm_q1, na.rm = TRUE),
      fwhm_q3 = stats::median(sub$fwhm_q3, na.rm = TRUE)
    )
  }))
  pooled_regions <- pooled_regions[order(pooled_regions$region), ]
  rownames(pooled_regions) <- NULL
  list(noise_sd = noise_sd, csnr = csnr, regions = pooled_regions)
}

## Region-aware counterpart to generate_synthetic_pool_1d(): pooled_stats
## here comes from pool_signal_stats_regions_1d() instead of
## pool_signal_stats_1d(), so each synthetic draw reproduces the real
## signals' own region-to-region density/fwhm variation via
## gen_synthetic_1d_regions() rather than one blended, uniform density.
#' @noRd
generate_synthetic_pool_1d_regions <- function(pooled_stats, y_list, n_synthetic = 10, peak_shape = "gaussian") {
  n <- round(stats::median(vapply(y_list, length, integer(1))))
  lapply(seq_len(n_synthetic), function(seed) {
    gen_synthetic_1d_regions(n = n, region_profile = pooled_stats$regions,
                              csnr = pooled_stats$csnr, peak_shape = peak_shape, seed = seed)
  })
}

## Greedily merges ADJACENT rows of a pool_signal_stats_regions_1d() regions
## table (left-to-right) so each merged group has an ESTIMATED peak count
## (density * span, in n_total's own point scale) of at least min_peaks --
## the region-level counterpart of tune_psalsa_params_1d()'s own n_total >=
## 15 threshold for a well-defined peak-area objective (see
## compute_peak_size_breaks_1d()). Addresses the root cause of a sparse
## region's tuning collapsing to a degenerate optimum (see
## tune_psalsa_region_params_1d()'s own comments): widening the region until
## it actually contains enough real peaks to constrain the search gives the
## objective genuine signal to work with, rather than relying on a strong
## prior to paper over too little data.
##
## A trailing group that still falls short after reaching the end is merged
## backward into the previous group rather than left under-informed. Only
## ADJACENT regions are merged (never skipping over a gap to reach a
## non-adjacent region with more peaks) -- growing a genuinely contiguous
## stretch of the spectrum is what "more data" means physically; splicing
## together two distant, unrelated regions would misrepresent the merged
## region's own span as uniformly dense when it isn't.
##
## Density/fwhm for a merged group are recombined as span-weighted mean
## density and min/max fwhm across its member rows (excluding empty rows,
## na.rm = TRUE) -- an approximation, but a representative one for building
## a synthetic pool over the merged span. A merged group with literally zero
## peaks across every member row (e.g. min_peaks can't be met anywhere
## because the whole tail of the spectrum is empty) keeps density 0/fwhm NA,
## same as an ordinary empty region -- still skipped by
## tune_psalsa_region_params_1d(), not force-tuned on nothing.
#' @noRd
merge_sparse_regions_1d <- function(regions, n_total, min_peaks = 15) {
  spans <- (regions$frac_hi - regions$frac_lo) * n_total
  est_peaks <- ifelse(is.na(regions$density), 0, regions$density * spans)

  groups <- integer(nrow(regions))
  g <- 1L
  acc <- 0
  for (i in seq_len(nrow(regions))) {
    groups[i] <- g
    acc <- acc + est_peaks[i]
    if (acc >= min_peaks && i < nrow(regions)) {
      g <- g + 1L
      acc <- 0
    }
  }
  ## Only merge a trailing group backward if it has SOME real signal but not
  ## enough (acc > 0): a genuinely empty trailing group (acc == 0) gains
  ## nothing from merging and would only dilute an already well-informed
  ## previous group's density and needlessly widen its span into empty
  ## territory -- left as its own (skippable, density-0) group instead.
  if (acc > 0 && acc < min_peaks && g > 1L) {
    groups[groups == g] <- g - 1L
  }

  merged_rows <- lapply(split(seq_len(nrow(regions)), groups), function(idx) {
    idx <- sort(idx)
    total_span <- sum(spans[idx])
    total_peaks <- sum(est_peaks[idx])
    has_peaks <- total_peaks > 0
    data.frame(
      frac_lo = regions$frac_lo[idx[1]], frac_hi = regions$frac_hi[idx[length(idx)]],
      density = if (has_peaks) total_peaks / total_span else 0,
      fwhm_q1 = if (has_peaks) min(regions$fwhm_q1[idx], na.rm = TRUE) else NA_real_,
      fwhm_q3 = if (has_peaks) max(regions$fwhm_q3[idx], na.rm = TRUE) else NA_real_
    )
  })
  merged <- do.call(rbind, merged_rows)
  merged <- merged[order(merged$frac_lo), ]
  merged$region <- seq_len(nrow(merged))
  rownames(merged) <- NULL
  merged[, c("region", "frac_lo", "frac_hi", "density", "fwhm_q1", "fwhm_q3")]
}

## For EACH region in pooled_regions$regions (pool_signal_stats_regions_1d()
## output), builds a synthetic pool of THAT region ALONE (length matched to
## its own share of n_total) and tunes lambda/p/k_mult against it
## independently via tune_psalsa_params_1d() -- one optimal parameter set
## per region, rather than one signal-wide compromise. A region with no
## observed peaks (density <= 0 or NA/degenerate fwhm -- see
## pool_signal_stats_regions_1d()) is skipped entirely: there's no local
## signal to tune against, and spatial_profile_1d() fills the gap smoothly
## from its neighbours' knots instead of guessing.
##
## theta0 anchors EACH region's own tune_psalsa_params_1d() ridge
## regularization (see that function's own comments for why the
## regularization exists at all: unconstrained Nelder-Mead can drift to a
## degenerate "disable peak protection" optimum on a single finite, noisy
## pool draw). tune_psalsa_params_1d()'s own default theta0 is a FIXED
## literature prior -- reasonable for a whole spectrum with plenty of peaks
## to outweigh it, but a single region typically has far fewer peaks, so its
## objective is flatter and more easily dominated by the SAME regularization
## pulling toward a prior that may not even suit this particular spectrum.
## Passing a theta0 built from THIS spectrum's own signal-wide tune (see
## tune_psalsa_spatial(), which computes it once and shares it across every
## region) anchors each region to what's typical for the spectrum at hand --
## "the pack" -- rather than a generic constant, and (still via the exact
## same regularization mechanism) keeps any one region's knot from swinging
## arbitrarily far from it just because that region's own peaks are too
## sparse to constrain the search on their own. NULL keeps
## tune_psalsa_params_1d()'s own literature-prior default (its previous,
## unregularized-by-context behaviour).
##
## region_lambda_k_weight/region_p_weight are the SAME regularization
## weights tune_psalsa_params_1d() uses for its own (whole-spectrum)
## lambda_k_weight/p_weight, but much stronger by default -- verified
## necessary, not just cautious: a genuinely sparse region (as few as 3
## peaks in its own synthetic pool) still collapsed to the SAME
## "disable peak protection" escape hatch (lambda->tiny, k_mult->huge,
## p->near p_max) even WITH a spectrum-specific theta0 anchor in place,
## at the global search's own default weights (0.02/0.4) -- confirmed
## reproducibly on a real region from the MTBLS242 dataset (region 20/20,
## 3 peaks/draw): lambda collapsed to ~400 vs. an anchor of ~2.8e7, and
## the resulting position-varying baseline interpolated the raw signal
## almost exactly in that region (baseline-corrected threshold collapsing
## to ~0 for every sample there, not just the ones the tuning was meant to
## fix). A region has far less data than the whole spectrum to resist that
## escape hatch, so the SAME ridge weight that suffices globally is
## comparably weaker locally; lambda_k_weight = 2 and p_weight = 10 (~100x
## and ~25x the global defaults) were verified to keep this same sparse
## region's tuned values close to its anchor (lambda within ~2x, k_mult
## within ~1.5x) across 5 independent synthetic pool draws.
##
## min_peaks: before tuning, ADJACENT regions with too few estimated peaks
## are merged (see merge_sparse_regions_1d()) so each region tuned below has
## at least this many real peaks to constrain its own search -- addressing
## the sparse-region collapse at its root cause (too little data) rather
## than relying only on the regularization above to paper over it. The two
## mechanisms are complementary, not redundant: merging can still leave a
## genuinely peak-poor stretch of spectrum under-informed (e.g. the very
## last group, or a spectrum with too few peaks anywhere to reach
## min_peaks), and the regularization above is what tune_psalsa_params_1d()
## itself already relies on even for a well-populated whole-spectrum search
## (see its own comments) -- it stays as a safety net regardless of how the
## region was sized.
##
## Returns a data frame (region, frac_mid, lambda, p, k_mult), one row per
## region that had enough signal to tune, or NULL if none did.
#' @noRd
tune_psalsa_region_params_1d <- function(pooled_regions, n_total, peak_shape = "gaussian",
                                          n_synthetic = 10, optim_maxit = 150, optim_reltol = 1e-6,
                                          theta0 = NULL, p_max = 0.05,
                                          region_lambda_k_weight = 2, region_p_weight = 10,
                                          min_peaks = 15) {
  regions <- merge_sparse_regions_1d(pooled_regions$regions, n_total = n_total, min_peaks = min_peaks)
  ## p_max is passed alongside theta0 (not on its own) because theta0's
  ## logit-transformed p component was encoded using THIS p_max -- decoding
  ## it with a different p_max inside tune_psalsa_params_1d() would silently
  ## misinterpret the anchor.
  extra_args <- if (is.null(theta0)) {
    list()
  } else {
    list(
      theta0 = theta0, p_max = p_max,
      lambda_k_weight = region_lambda_k_weight, p_weight = region_p_weight
    )
  }
  rows <- lapply(seq_len(nrow(regions)), function(r) {
    reg <- regions[r, ]
    if (is.na(reg$density) || reg$density <= 0 ||
      is.na(reg$fwhm_q1) || is.na(reg$fwhm_q3) || reg$fwhm_q3 <= 0) {
      return(NULL)
    }
    reg_n <- max(10, round((reg$frac_hi - reg$frac_lo) * n_total))
    pool_r <- lapply(seq_len(n_synthetic), function(seed) {
      gen_synthetic_1d(n = reg_n, density = reg$density, fwhm_range = c(reg$fwhm_q1, reg$fwhm_q3),
                        csnr = pooled_regions$csnr, peak_shape = peak_shape, seed = seed)
    })
    tuned_r <- do.call(tune_psalsa_params_1d, c(
      list(pool_r, optim_maxit = optim_maxit, optim_reltol = optim_reltol), extra_args
    ))
    data.frame(
      region = reg$region, frac_mid = (reg$frac_lo + reg$frac_hi) / 2,
      lambda = tuned_r$lambda, p = tuned_r$p, k_mult = tuned_r$k_mult
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

## n knot positions (frac_mid * n) + knot values -> a smooth length-n
## position profile, interpolated in a TRANSFORMED space chosen so the
## inverse transform always lands back in the parameter's valid range no
## matter how the spline over/undershoots between knots:
##   "log"      (lambda, k_mult: positive, no upper bound) -- exp() of
##              anything real is always > 0.
##   "logit"    (p: bounded in (0, p_max)) -- plogis() of anything real
##              always lands in (0, 1), scaled to (0, p_max).
## A natural cubic spline is used (not linear) so the resulting lambda/p/k
## profile is itself smooth -- a piecewise-linear profile would put a kink
## in the penalty at every region boundary, which is exactly the
## discontinuity this design is meant to avoid.
#' @noRd
spatial_profile_1d <- function(n, frac_mid, values, transform = c("log", "logit", "identity"), p_max = 0.05) {
  transform <- match.arg(transform)
  y_knots <- switch(transform,
    log = log(values),
    logit = stats::qlogis(values / p_max),
    identity = values
  )
  z <- if (length(frac_mid) < 2) {
    rep(y_knots[1], n)
  } else {
    stats::spline(x = frac_mid * n, y = y_knots, xout = seq_len(n), method = "natural")$y
  }
  switch(transform,
    log = exp(z),
    logit = stats::plogis(z) * p_max,
    identity = z
  )
}

## Stretches/shrinks a length-n_src position profile onto n_target points by
## linear interpolation over the shared FRACTIONAL axis (0..1) -- used when
## applying one tuned spatial profile (built at the pooled median length) to
## an individual spectrum of a different length. rule = 2 (nearest-edge
## extrapolation) only matters at the very ends, if n_target's fractional
## grid extends fractionally past n_src's.
#' @noRd
resample_profile_1d <- function(profile, n_target) {
  n_src <- length(profile)
  if (n_src == n_target) return(profile)
  stats::approx(
    x = (seq_len(n_src) - 0.5) / n_src, y = profile,
    xout = (seq_len(n_target) - 0.5) / n_target, rule = 2
  )$y
}

## Nelder-Mead-tunes psalsa(lambda, p, k) against a synthetic pool with known
## ground truth. Balances recovery across small/medium/large true-area tiers
## (each capped at 1000% error, tier dropped if fewer than 2 peaks) rather
## than optimizing a single pooled error, so a few large peaks can't dominate
## the objective at the expense of small ones. Falls back to floor_rmse_1d()
## (still against the pool's own known ground truth) when the pool has too
## few ISOLATED peaks to define three size tiers meaningfully.
##
## p is capped well below 0.5 (p_max) and ridge-penalized toward the
## literature-default prior (theta0), rather than left to range freely in
## (0, 0.5) -- found (via repeated, reproducible failures on very ordinary
## single-peak test signals, both with the peak-area objective active and
## with plenty of isolated peaks available) that unconstrained Nelder-Mead
## can still drift to p~0.5, disabling PSALSA's own asymmetric peak
## protection, on a single finite/noisy synthetic pool draw where that
## region happens to look locally attractive -- even though nothing forces
## it there and the correct objective is in play. p=0.5 isn't just an
## extreme choice, it's the exact point where p==1-p and PSALSA's asymmetry
## vanishes mathematically, so foreclosing that region (not merely
## discouraging it) is justified on first principles: literature/typical p
## values are always small (~0.001-0.05) anyway. The ridge penalty makes
## ANY large excursion from the sane prior costly unless the data
## genuinely earns it.
##
## Regularizing p ALONE isn't enough: verified that once p is pinned near
## its prior, the optimizer just relocates the exact same "disable
## asymmetric protection" outcome to a DIFFERENT escape hatch instead --
## lambda collapsing toward 0 together with k exploding (both, independent
## of p, also eliminate any real peak resistance). Peak protection can be
## disabled via p->0.5, OR via lambda->0, OR via k->Inf, OR any mix, so all
## three get SOME ridge treatment, anchored at their own theta0 component --
## but NOT with the same weight. p's "right" value is genuinely universal
## (should almost always be small), so it gets a strong anchor
## (p_weight=0.4). lambda's "right" value is NOT universal -- it depends on
## THIS signal's own baseline curvature scale, which varies signal to
## signal -- so anchoring it as strongly as p systematically OVERSMOOTHS
## whenever the true scale differs from the prior (verified: a real signal
## whose baseline curvature genuinely needed a ~20x smaller lambda than the
## 1e7 prior got dragged most of the way back to it anyway, RMSE 3.2 vs an
## achievable ~0.5, even though nothing was numerically unstable). A much
## lighter shared weight (lambda_k_weight=0.02) on lambda/k_mult still
## reliably blocks the same joint runaway (verified across 8 random noise
## seeds: too light, e.g. 0.005, and lambda collapses toward 0 with the
## SAME peak-climbing failure as no regularization at all; 0.02 stays
## stable and sane in every seed while still letting lambda move ~20x off
## the prior when the data genuinely supports it).
#' @noRd
tune_psalsa_params_1d <- function(pool, p_max = 0.05, p_weight = 0.4, lambda_k_weight = 0.02,
                                   theta0 = c(log(1e7), stats::qlogis(0.001 / p_max), log(15)),
                                   optim_maxit = 150, optim_reltol = 1e-6) {
  size_breaks <- compute_peak_size_breaks_1d(pool)
  n_total <- sum(vapply(pool, function(sig) nrow(sig$peak_info), integer(1)))
  use_peak_area <- all(is.finite(size_breaks)) &&
    length(unique(c(-Inf, size_breaks, Inf))) == 4 && n_total >= 15

  regularization <- function(theta) {
    p_weight * (theta[2] - theta0[2])^2 + lambda_k_weight * ((theta[1] - theta0[1])^2 + (theta[3] - theta0[3])^2)
  }

  obj <- if (use_peak_area) {
    function(theta) {
      lambda <- exp(theta[1]); p <- stats::plogis(theta[2]) * p_max; k_mult <- exp(theta[3])
      pae_rows <- list()
      for (sig in pool) {
        z <- tryCatch(as.vector(psalsa(sig$y, lambda = lambda, p = p, k = k_mult * noise_est_1d(sig$y))$baseline),
                      error = function(e) NULL)
        if (is.null(z) || any(!is.finite(z))) next
        pae <- peak_area_errors_1d(sig$y, z, sig$peak_info, sig$peaks, sig$peak_shape)
        pae_rows[[length(pae_rows) + 1]] <- cbind(pae, size_tier = classify_peak_size_1d(pae$area, size_breaks))
      }
      if (!length(pae_rows)) return(1e6)
      all_pae <- do.call(rbind, pae_rows)
      tier_means <- sapply(c("small", "medium", "large"), function(tr) {
        sub <- all_pae[all_pae$size_tier == tr & !is.na(all_pae$size_tier), ]
        if (nrow(sub) < 2) return(NA_real_)
        mean(pmin(sub$abs_pct_error, 1000))
      })
      if (all(is.na(tier_means))) return(1e6)
      mean(tier_means, na.rm = TRUE) + regularization(theta)
    }
  } else {
    function(theta) {
      lambda <- exp(theta[1]); p <- stats::plogis(theta[2]) * p_max; k_mult <- exp(theta[3])
      errs <- vapply(pool, function(sig) {
        z <- tryCatch(as.vector(psalsa(sig$y, lambda = lambda, p = p, k = k_mult * noise_est_1d(sig$y))$baseline),
                      error = function(e) NULL)
        if (is.null(z) || any(!is.finite(z))) return(NA_real_)
        floor_rmse_1d(z, sig$baseline, sig$peaks, sig$sigma, mean(sig$baseline))
      }, numeric(1))
      if (all(is.na(errs))) return(1e6)
      mean(errs, na.rm = TRUE) + regularization(theta)
    }
  }

  opt <- stats::optim(theta0, obj, method = "Nelder-Mead", control = list(maxit = optim_maxit, reltol = optim_reltol))
  ## k_mult, not an absolute k: the pool's synthetic signals don't share the
  ## real input's noise scale, so the caller (tune_psalsa()) turns this into
  ## an absolute k using the REAL input's own noise_est_1d().
  list(lambda = exp(opt$par[1]), p = stats::plogis(opt$par[2]) * p_max, k_mult = exp(opt$par[3]))
}
