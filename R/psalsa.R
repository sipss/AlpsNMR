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
#' @param lambda Smoothing parameter (generally `1e5` - `1e8`).
#' @param p Asymmetry parameter.
#' @param k Peak height parameter, controlling how strongly the weights of
#'   points above the baseline decay with their height. Usually about 5% of the
#'   maximum intensity. When `k = -1` (the default) it is set to one twentieth
#'   of the maximum intensity of each spectrum.
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
  penalty <- lambda * diff2_penalty(length(y))
  psalsa_core(y, function(y, w) whit1d(y, penalty, w), p, k, maxit, k_epsilon)
}

#' 1D weighted Whittaker smoother
#'
#' Smooths a vector `y` with per-element weights `w` by solving the penalised
#' least-squares system `(diag(w) + penalty) s = w * y`, where `penalty` is the
#' pre-built second-order difference penalty (`lambda * t(D2) %*% D2`, see
#' [diff2_penalty()]) -- the classic second-order Whittaker smoother, solved
#' directly via sparse Cholesky. Used both by [psalsa()] directly and, along
#' the processing axis, by [psalsa2d()].
#'
#' @param y Numeric vector to smooth.
#' @param penalty Sparse penalty matrix (`lambda`-scaled) from [diff2_penalty()].
#' @param w Numeric vector of weights, same length as `y`.
#' @return A numeric vector with the smoothed signal.
#' @noRd
whit1d <- function(y, penalty, w) {
  a <- Matrix::forceSymmetric(Matrix::Diagonal(x = w) + penalty)
  as.vector(Matrix::solve(a, w * y))
}

#' PSALSA asymmetric reweighting core
#'
#' Shared iteration used by both the 1D ([psalsa()]) and 2D ([psalsa2d()])
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
  if (k == -1) {
    k <- max(y) / 20
  }
  k <- max(k, k_epsilon * diff(range(y)))

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
    w[d_geq] <- p * exp(-d[d_geq] / k)
    w[!d_geq] <- 1 - p

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
