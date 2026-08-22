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
psalsa <- function(spectra, lambda = 1e+07, p = 0.001, k = -1, maxit = 25, k_epsilon = 1e-6) {
  if (maxit < 1) {
    stop("maxit smaller than 1")
  }

  if (is.matrix(spectra)) {
    estbaseline <- 0 * spectra
    for (i in seq_len(nrow(spectra))) {
      estbaseline[i, ] <- psalsa_one(spectra[i, ], lambda, p, k, maxit, k_epsilon)
    }
  } else {
    estbaseline <- psalsa_one(spectra, lambda, p, k, maxit, k_epsilon)
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
psalsa_one <- function(y, lambda = 1e+07, p = 0.001, k = -1, maxit = 25, k_epsilon = 1e-6) {
  ## Scalar lambda keeps the original code path byte-for-byte (lambda * the
  ## unweighted penalty), rather than routing it through
  ## diff2_penalty_weighted()'s sqrt-then-crossprod construction, which is
  ## mathematically equivalent but not guaranteed bit-identical.
  penalty <- if (length(lambda) == 1) {
    lambda * diff2_penalty(length(y))
  } else {
    diff2_penalty_weighted(length(y), lambda)
  }
  psalsa_core(y, function(y, w) whit1d(y, penalty, w), p, k, maxit, k_epsilon)
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
psalsa_core <- function(y, smoother, p = 0.001, k = -1, maxit = 25, k_epsilon = 1e-6) {
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
    s <- smoother(y, w)

    d_geq_old <- d_geq
    d <- y - s
    d_geq <- d >= 0
    w[d_geq] <- p[d_geq] * exp(-d[d_geq] / k[d_geq])
    w[!d_geq] <- 1 - p[!d_geq]

    # converged once the set of points above the baseline no longer changes
    if (all(d_geq == d_geq_old)) {
      break
    }
  }

  s
}

diff2_penalty <- function(n) {
  d2 <- Matrix::bandSparse(
    n - 2, n,
    k = c(0, 1, 2),
    diagonals = list(rep(1, n - 2), rep(-2, n - 2), rep(1, n - 2))
  )
  Matrix::crossprod(d2)
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
