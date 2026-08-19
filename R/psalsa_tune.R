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
#' @param peak_shape Either `"gaussian"` or `"gex"` (an asymmetric,
#'   exponentially-modified peak shape typical of chromatography). Peak shape
#'   is often known from the instrument/technique used to acquire `y`, so it
#'   is left as an explicit argument rather than inferred from the data.
#' @param n_synthetic Number of synthetic spectra generated for the search.
#' @param optim_maxit Maximum number of [stats::optim()] (Nelder-Mead)
#'   iterations.
#' @param optim_reltol [stats::optim()] relative convergence tolerance.
#'
#' @family baseline removal functions
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
tune_psalsa <- function(y, peak_shape = c("gaussian", "gex"), n_synthetic = 10,
                         optim_maxit = 150, optim_reltol = 1e-6) {
  peak_shape <- match.arg(peak_shape)
  y_list <- if (is.list(y)) y else list(y)
  pooled_stats <- pool_signal_stats_1d(y_list)
  pool <- generate_synthetic_pool_1d(pooled_stats, y_list, n_synthetic = n_synthetic, peak_shape = peak_shape)
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

## Peak-area recovery: for each peak in peak_info, integrates the baseline-
## corrected signal (y - z_est) over that peak's own window and compares to
## its true area.
#' @noRd
peak_area_errors_1d <- function(y, z_est, peak_info) {
  if (is.null(peak_info) || nrow(peak_info) == 0) return(peak_info)
  corrected <- as.vector(y) - as.vector(z_est)
  est_area <- vapply(seq_len(nrow(peak_info)), function(i) {
    sum(corrected[peak_info$lo[i]:peak_info$hi[i]])
  }, numeric(1))
  pct_err <- 100 * (est_area - peak_info$area) / peak_info$area
  cbind(peak_info, est_area = est_area, pct_error = pct_err, abs_pct_error = abs(pct_err),
        isolated = peak_is_isolated_1d(peak_info))
}

## Small/medium/large tiers by true peak area, using cut points computed once
## from the pool's own isolated-peak population, so tiers are comparable
## across the different synthetic signals scored by the same tuning run.
#' @noRd
compute_peak_size_breaks_1d <- function(signals) {
  areas <- unlist(lapply(signals, function(s) {
    pi <- s$peak_info
    if (is.null(pi) || nrow(pi) == 0) return(NULL)
    pi$area[peak_is_isolated_1d(pi)]
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
## objective when the pool has too few isolated peaks for peak-area scoring.
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
  n_iso <- sum(vapply(pool, function(sig) sum(peak_is_isolated_1d(sig$peak_info)), integer(1)))
  use_peak_area <- all(is.finite(size_breaks)) &&
    length(unique(c(-Inf, size_breaks, Inf))) == 4 && n_iso >= 15

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
        pae <- peak_area_errors_1d(sig$y, z, sig$peak_info)
        pae_rows[[length(pae_rows) + 1]] <- cbind(pae, size_tier = classify_peak_size_1d(pae$area, size_breaks))
      }
      if (!length(pae_rows)) return(1e6)
      iso <- do.call(rbind, pae_rows); iso <- iso[iso$isolated, ]
      tier_means <- sapply(c("small", "medium", "large"), function(tr) {
        sub <- iso[iso$size_tier == tr & !is.na(iso$size_tier), ]
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
