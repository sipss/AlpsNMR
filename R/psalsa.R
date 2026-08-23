#' Baseline estimation with Peaked Signal's Asymmetric Least Squares
#'
#' Estimates the baseline of one or more spectra using the PSALSA
#' (Peaked Signal's Asymmetric Least Squares Algorithm) method. PSALSA is a
#' variant of asymmetric least squares that limits the influence of intense
#' peaks on the estimated baseline, which improves the result for signals
#' containing sharp, high peaks.
#'
#' Like ordinary asymmetric least squares, PSALSA assigns different weights to
#' the points above and below an iteratively estimated baseline: the asymmetry
#' parameter `p` (`0 <= p <= 1`) is the weight for points below the baseline,
#' whereas points above it normally receive weight `1 - p`. The difference is
#' that for points above the baseline the weight decays exponentially with the
#' height of the point above the current estimate, controlled by `k`. This
#' prevents intense peaks from pulling the baseline upwards. The parameter
#' `lambda` controls the amount of smoothing: the larger it is, the smoother
#' the baseline will be. Iteration stops once the set of points above the
#' baseline no longer changes, or after `maxit` iterations.
#'
#' The baseline is fitted with a weighted second-order Whittaker smoother.
#'
#' @param spectra Either a numeric vector containing one spectrum, or a matrix
#'   containing one spectrum per row.
#' @param lambda Smoothing parameter (generally `1e5` - `1e8`). Either a single
#'   value (used for the whole spectrum) or a numeric vector the same length
#'   as one spectrum, giving a position-varying smoothing profile -- see
#'   `tune_psalsa_spatial()`, which builds such a profile from per-region
#'   tuning. When `spectra` is a matrix, the same profile is used for every
#'   row.
#' @param p Asymmetry parameter. Either a single value or a position-varying
#'   vector, as for `lambda`.
#' @param k Peak height parameter, controlling how strongly the weights of
#'   points above the baseline decay with their height. Usually about 5% of the
#'   maximum intensity. When `k = -1` (the default) it is set to one twentieth
#'   of the maximum intensity of each spectrum. Either a single value or a
#'   position-varying vector, as for `lambda` (the `k = -1` auto-default only
#'   applies when `k` is a single value).
#' @param maxit Maximum number of iterations.
#' @param k_epsilon Floors `k` at `k_epsilon * diff(range(spectra))`, so an
#'   estimated or supplied `k` that is (near) zero can't make the asymmetric
#'   weight `exp(-d/k)` collapse to 0 for essentially any positive residual --
#'   which happens for a genuinely near-noise-free signal, since `k` is
#'   usually derived from an estimated noise level (see [tune_psalsa()]),
#'   silently degenerating the baseline toward the signal's minimum. Scaled by
#'   the data's own range (rather than a fixed absolute value) so the floor
#'   stays proportionate regardless of the signal's amplitude.
#' @param damping Under-relaxation factor in `(0, 1]` for the asymmetric
#'   weight update: each iteration's new weight is `damping * w_target + (1 -
#'   damping) * w_old` instead of replacing `w` outright. `damping = 1` (the
#'   default) is the original, undamped update. Because the reweighting
#'   depends on a hard threshold (whether each point sits above or below the
#'   current baseline estimate), a full-step update can occasionally
#'   overshoot near a flat/quiet region -- many points sit close to that
#'   threshold there, so a large-enough step in the fitted curve flips a big
#'   batch of them at once, which then swings the next fitted curve enough to
#'   flip a comparable batch back, an oscillation that can persist right up
#'   to `maxit` instead of settling. `damping < 1` shrinks how much the fitted
#'   curve moves per iteration, directly shrinking how many points can cross
#'   the threshold in one step. Verified on real data to eliminate this
#'   oscillation (see `psalsa_core()`'s own comments for the failure mode);
#'   convergence takes correspondingly more iterations, so `maxit` may need
#'   raising alongside a `damping` below 1.
#' @param weight_method `"threshold"` (the default) is the original weight
#'   update (`p * exp(-d/k)` above the current baseline, `1 - p` below it,
#'   split by the hard `d >= 0` threshold). `"smooth"` instead uses
#'   `compute_psalsa_weights()`: a piecewise, C1-continuous weight curve with
#'   a central noise band `[-s, s]` (cubic Hermite-interpolated between the
#'   two boundary weights, matching the exponential pieces' own slopes there)
#'   instead of a hard threshold at `d = 0`. Two effects, addressing two
#'   distinct issues found with `"threshold"`: (1) no discontinuity at `d =
#'   0` to oscillate across (may reduce or remove the need for `damping <
#'   1`); (2) points near the middle of the noise band get a weight between
#'   `p` and `w_max` rather than `p` itself, avoiding the systematic
#'   downward bias `"threshold"` introduces even in a genuinely flat/quiet
#'   region (there, `p` alone -- not `k` -- sets a floor on how low a
#'   residual's weight can go, however small the residual, since `d >= 0`
#'   alone triggers the `p`-branch). Only supports a single (scalar) `p`/`k`
#'   for now, not `tune_psalsa_spatial()`'s position-varying profiles.
#' @param s,c_noise,w_max Only used when `weight_method = "smooth"` -- see
#'   `compute_psalsa_weights()`. `s` should be sized to the signal's own
#'   noise amplitude (e.g. a small multiple of its estimated noise standard
#'   deviation), not left at its exploratory default of `500` for a signal
#'   whose real noise level is very different.
#'
#' @return A list with two elements, each with the same dimensions as
#'   `spectra`:
#'   \describe{
#'     \item{`baseline`}{the estimated baseline.}
#'     \item{`corrected`}{the baseline-corrected signal, i.e.
#'       `spectra - baseline`.}
#'   }
#'
#' @references
#' Oller-Moreno, S., Pardo, A., Jimenez-Soto, J. M., Samitier, J., Marco, S.
#' (2014). "Adaptive Asymmetric Least Squares baseline estimation for analytical
#' instruments". 2014 IEEE 11th International Multi-Conference on Systems,
#' Signals & Devices (SSD14), 1-5. \doi{10.1109/SSD.2014.6808837}
#'
#' @examples
#' x <- seq_len(200)
#' baseline <- 10 + 5 * sin(x / 40)
#' peak <- 80 * exp(-((x - 120)^2) / (2 * 4^2))
#' y <- baseline + peak
#'
#' result <- AlpsNMR:::psalsa(y)
#' plot(y, type = "l")
#' lines(result$baseline, col = "red")
#'
psalsa <- function(spectra, lambda = 1e+07, p = 0.001, k = -1, maxit = 25, k_epsilon = 1e-6, damping = 1,
                    weight_method = c("threshold", "smooth"), s = 500, c_noise = 0.1, w_max = 0.9) {
  weight_method <- match.arg(weight_method)
  if (maxit < 1) {
    stop("maxit smaller than 1")
  }

  if (is.matrix(spectra)) {
    estbaseline <- 0 * spectra
    for (i in seq_len(nrow(spectra))) {
      estbaseline[i, ] <- psalsa_one(spectra[i, ], lambda, p, k, maxit, k_epsilon, damping, weight_method, s, c_noise, w_max)
    }
  } else {
    estbaseline <- psalsa_one(spectra, lambda, p, k, maxit, k_epsilon, damping, weight_method, s, c_noise, w_max)
  }

  list(baseline = estbaseline, corrected = spectra - estbaseline)
}

#' PSALSA baseline for a single spectrum
#'
#' Internal worker used by [psalsa()] to estimate the baseline of one spectrum.
#'
#' @param y Numeric vector with one spectrum.
#' @inheritParams psalsa
#' @return A numeric vector with the estimated baseline.
#' @noRd
psalsa_one <- function(y, lambda = 1e+07, p = 0.001, k = -1, maxit = 25, k_epsilon = 1e-6, damping = 1,
                        weight_method = "threshold", s = 500, c_noise = 0.1, w_max = 0.9) {
  ## Scalar lambda keeps the original code path byte-for-byte (lambda * the
  ## unweighted penalty), rather than routing it through
  ## diff2_penalty_weighted()'s sqrt-then-crossprod construction, which is
  ## mathematically equivalent but not guaranteed bit-identical.
  penalty <- if (length(lambda) == 1) {
    lambda * diff2_penalty(length(y))
  } else {
    diff2_penalty_weighted(length(y), lambda)
  }
  psalsa_core(y, function(y, w) whit1d(y, penalty, w), p, k, maxit, k_epsilon, damping, weight_method, s, c_noise, w_max)
}

#' 1D weighted Whittaker smoother
#'
#' Smooths a vector `y` with per-element weights `w` by solving the penalised
#' least-squares system `(diag(w) + penalty) s = w * y`, where `penalty` is the
#' pre-built second-order difference penalty (`lambda * t(D2) %*% D2`, see
#' `diff2_penalty()`) -- the classic second-order Whittaker smoother, solved
#' directly via sparse Cholesky. Used both by [psalsa()] directly and, along
#' the processing axis, by `psalsa2d()`.
#'
#' @param y Numeric vector to smooth.
#' @param penalty Sparse penalty matrix (`lambda`-scaled) from `diff2_penalty()`.
#' @param w Numeric vector of weights, same length as `y`.
#' @return A numeric vector with the smoothed signal.
#' @noRd
whit1d <- function(y, penalty, w) {
  a <- Matrix::forceSymmetric(Matrix::Diagonal(x = w) + penalty)
  as.vector(Matrix::solve(a, w * y))
}

#' PSALSA asymmetric reweighting core
#'
#' Shared iteration used by both the 1D ([psalsa()]) and 2D (`psalsa2d()`)
#' baseline estimators. The only difference between the two is the smoother, so
#' it is passed in as a function `smoother(y, w)` that returns the smoothed
#' signal (of the same shape as `y`) for the weight array `w`. Because every
#' operation in the loop is elementwise, the same code works for a numeric
#' vector (1D) or a matrix (2D image).
#'
#' @param y Numeric vector or matrix with the signal.
#' @param smoother Function of `(y, w)` returning the smoothed signal.
#' @inheritParams psalsa
#' @return The estimated baseline, with the same shape as `y`.
#' @noRd
psalsa_core <- function(y, smoother, p = 0.001, k = -1, maxit = 25, k_epsilon = 1e-6, damping = 1,
                         weight_method = "threshold", s = 500, c_noise = 0.1, w_max = 0.9) {
  if (weight_method == "smooth" && (length(p) != 1 || length(k) != 1)) {
    stop("weight_method = \"smooth\" only supports a single (scalar) p/k, not a position-varying profile.")
  }

  ## The k = -1 auto-default only makes sense for a single value -- a
  ## position-varying k must be supplied explicitly by the caller (e.g.
  ## tune_psalsa_spatial()'s k_mult_profile * noise_sd).
  if (length(k) == 1 && k == -1) {
    k <- max(y) / 20
  }
  k <- pmax(k, k_epsilon * diff(range(y)))

  ## Broadcast p/k to the signal's own length up front, so the reweighting
  ## step below can safely subset them alongside d[d_geq]/d[!d_geq] whether
  ## they were passed as a single value (recycled, byte-identical to the old
  ## scalar-only code) or as a position-varying vector -- subsetting a
  ## not-yet-broadcast vector with a longer logical mask would misalign it.
  p <- rep_len(p, length(y))
  k <- rep_len(k, length(y))

  # weights, same shape as the signal, all ones so the first pass is unweighted
  w <- y
  w[] <- 1
  # sign pattern of the residuals; all FALSE so the first iteration always runs
  d_geq <- (w < 0)

  for (it in seq_len(maxit)) {
    baseline <- smoother(y, w)

    d_geq_old <- d_geq
    d <- y - baseline
    d_geq <- d >= 0
    w_target <- if (weight_method == "smooth") {
      compute_psalsa_weights(d, s = s, p = p[1], c_noise = c_noise, k = k[1], w_max = w_max)
    } else {
      w_target <- w
      w_target[d_geq] <- p[d_geq] * exp(-d[d_geq] / k[d_geq])
      w_target[!d_geq] <- 1 - p[!d_geq]
      w_target
    }
    ## damping = 1 (the default) makes this an exact replacement, identical
    ## to the original undamped update -- see the `damping` parameter's own
    ## docs in psalsa() for why a smaller value is sometimes needed.
    w <- if (damping == 1) w_target else damping * w_target + (1 - damping) * w

    # converged once the set of points above the baseline no longer changes
    if (all(d_geq == d_geq_old)) {
      break
    }
  }

  baseline
}

diff2_penalty <- function(n) {
  d2 <- Matrix::bandSparse(
    n - 2, n,
    k = c(0, 1, 2),
    diagonals = list(rep(1, n - 2), rep(-2, n - 2), rep(1, n - 2))
  )
  Matrix::crossprod(d2)
}

#' Asymmetric Weight Function (C1 Continuous) for psalsa
#'
#' Calculates the weights w(d) for the residuals (y - z) using a continuously
#' differentiable piecewise formulation, based on cubic Hermite interpolation
#' for the noise region.
#'
#' @param d Numeric vector of residuals (signal y - baseline z).
#' @param s Noise threshold (half the width of the central band).
#' @param p Weight assigned at the boundaries of the noise band (|d| = s).
#' @param c_noise Approximate minimum weight for the desired noise (helps define mL).
#' @param k Tau constant (decay scale) for peaks (d > s).
#' @param w_max Asymptotic maximum penalty for severe underestimates (d < -s).
#'
#' @return A numeric vector of weights of the same length as `d`.
#' @export
compute_psalsa_weights <- function(d, s = 500, p = 0.3, c_noise = 0.1, k = 50, w_max = 0.9) {
  ## The left-region (d < -s) exponential is only guaranteed to saturate
  ## toward w_max (rather than diverge) when gamma_L > 0, which requires
  ## m_L < 0, which requires p > c_noise -- verified directly: with c_noise
  ## > p (e.g. PSALSA's typical p, often << 0.1), gamma_L flips negative and
  ## the exponential grows unboundedly for any point sitting meaningfully
  ## below the current baseline estimate, reaching -Inf within a few
  ## reweighting iterations on real data.
  if (c_noise > p) {
    stop("compute_psalsa_weights(): c_noise must be <= p, or the left-region weight diverges instead of saturating to w_max.")
  }

  # 1. Define slopes at the boundaries
  # m_L: Left slope at d = -s (approximated based on a virtual parabola)
  # This allows the user to avoid guessing m_L directly.
  m_L <- -2 * (p - c_noise) / s

  # m_R: Right slope at d = s (dictated by the derivative of the exponential)
  m_R <- -p / k

  # Pre-allocate results vector
  w <- numeric(length(d))

  # Indices for each region
  idx_left  <- d < -s
  idx_right <- d > s
  idx_mid   <- !idx_left & !idx_right

  # 2. Left Region (d < -s): High asymptotic penalty
  if (any(idx_left)) {
    # Exponent constant: -m_L / (w_max - p)
    # Ensures that the derivative evaluated at d = -s is exactly m_L
    gamma_L <- -m_L / (w_max - p)
    w[idx_left] <- w_max - (w_max - p) * exp(gamma_L * (d[idx_left] + s))
  }

  # 3. Right Region (d > s): Exponential decay to ignore peaks
  if (any(idx_right)) {
    w[idx_right] <- p * exp(-(d[idx_right] - s) / k)
  }

  # 4. Central Region (-s <= d <= s): Cubic Hermite Polynomial
  if (any(idx_mid)) {
    # Normalize d to the interval [0, 1]
    t <- (d[idx_mid] + s) / (2 * s)

    # Pre-compute powers
    t2 <- t^2
    t3 <- t^3

    # Hermite basis functions
    h00 <-  2 * t3 - 3 * t2 + 1
    h10 <-      t3 - 2 * t2 + t
    h01 <- -2 * t3 + 3 * t2
    h11 <-      t3 -     t2

    # The interval width in 'd' is (2*s).
    # By the chain rule (dt/dd = 1/(2s)), we must scale
    # the tangent vectors m_L and m_R by multiplying them by (2*s).
    scale_factor <- 2 * s

    # Assemble the Hermite polynomial
    w[idx_mid] <- p * h00 +
                 (m_L * scale_factor) * h10 +
                 p * h01 +
                 (m_R * scale_factor) * h11
  }

  return(w)
}

## Position-varying counterpart to `lambda * diff2_penalty(n)`: builds
## `t(D2) %*% diag(lambda_rows) %*% D2` instead of a single scalar-weighted
## penalty, so the smoothness constraint can be locally loosened or
## tightened rather than applying one constant lambda everywhere. Same
## pentadiagonal sparsity pattern as the scalar case regardless (each row of
## D2 only touches 3 neighbouring points, so diag(lambda_rows) never adds
## fill-in) -- no memory blowup from going position-varying.
##
## lambda may be length n (a full position profile, e.g. from
## spatial_profile_1d() -- the natural convention shared with p/k, which are
## used at every point) or length n - 2 (one weight per D2 row directly), or
## a scalar (recycled). A length-n profile is sliced to its interior
## 2:(n - 1) points, matching each D2 row's own center point.
##
## crossprod(diag(sqrt(w)) %*% d2) == t(d2) %*% diag(w) %*% d2: the sqrt
## split is what lets Matrix::crossprod() build the weighted penalty
## directly (crossprod computes t(A) %*% A), without ever materializing a
## dense diag(w) or losing sparsity.
#' @noRd
diff2_penalty_weighted <- function(n, lambda) {
  d2 <- Matrix::bandSparse(
    n - 2, n,
    k = c(0, 1, 2),
    diagonals = list(rep(1, n - 2), rep(-2, n - 2), rep(1, n - 2))
  )
  row_w <- if (length(lambda) == n) lambda[2:(n - 1)] else rep_len(lambda, n - 2)
  Matrix::crossprod(Matrix::Diagonal(x = sqrt(row_w)) %*% d2)
}
