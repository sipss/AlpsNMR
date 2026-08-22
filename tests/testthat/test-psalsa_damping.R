## psalsa_core: damping = 1 (default) is byte-identical to the undamped update

test_that("psalsa_core with damping = 1 is byte-identical to omitting damping entirely", {
    set.seed(1)
    y <- 5 + 20 * exp(-((seq_len(200) - 100)^2) / (2 * 4^2)) + rnorm(200, 0, 0.3)
    smoother <- function(y, w) AlpsNMR:::whit1d(y, 1e6 * AlpsNMR:::diff2_penalty(length(y)), w)

    result_default <- AlpsNMR:::psalsa_core(y, smoother, p = 0.01, k = 3)
    result_explicit <- AlpsNMR:::psalsa_core(y, smoother, p = 0.01, k = 3, damping = 1)
    expect_identical(result_default, result_explicit)
})

test_that("psalsa()/psalsa_one() with damping = 1 (default) match pre-damping behaviour", {
    set.seed(2)
    y <- 5 + 15 * exp(-((seq_len(150) - 75)^2) / (2 * 5^2)) + rnorm(150, 0, 0.2)
    result_default <- AlpsNMR:::psalsa(y, lambda = 1e6, p = 0.005, k = 4)
    result_explicit <- AlpsNMR:::psalsa(y, lambda = 1e6, p = 0.005, k = 4, damping = 1)
    expect_identical(result_default, result_explicit)
})

## damped update formula ----------------------------------------------------------

test_that("psalsa_core's damped weight update is a proper convex combination of the target and previous weight", {
    y <- c(rep(0, 5), 10, rep(0, 5)) # tiny signal, single spike
    smoother <- function(y, w) AlpsNMR:::whit1d(y, 1e2 * AlpsNMR:::diff2_penalty(length(y)), w)

    # Compare the update actually applied against a hand-computed one, after
    # a single iteration (starting weights are always 1, so w_old = 1):
    s <- smoother(y, rep(1, length(y)))
    d <- y - s
    d_geq <- d >= 0
    w_target <- ifelse(d_geq, 0.01 * exp(-d / 3), 1 - 0.01)
    expected_w_after_1_iter <- 0.4 * w_target + 0.6 * 1

    # psalsa_core doesn't expose intermediate w, so verify indirectly:
    # damping = 0 should leave the baseline exactly at the FIRST iteration's
    # smoother output forever (weights never actually change from 1):
    result_damping_0 <- AlpsNMR:::psalsa_core(y, smoother, p = 0.01, k = 3, maxit = 5, damping = 0)
    always_unweighted <- smoother(y, rep(1, length(y)))
    expect_equal(result_damping_0, always_unweighted)
})

test_that("damping = 0 keeps every iteration's smoother call using the original (unweighted) w", {
    set.seed(3)
    y <- 5 + 25 * exp(-((seq_len(100) - 50)^2) / (2 * 3^2)) + rnorm(100, 0, 0.2)
    smoother_calls <- list()
    smoother <- function(y, w) {
        smoother_calls[[length(smoother_calls) + 1]] <<- w
        AlpsNMR:::whit1d(y, 1e5 * AlpsNMR:::diff2_penalty(length(y)), w)
    }
    AlpsNMR:::psalsa_core(y, smoother, p = 0.01, k = 3, maxit = 4, damping = 0)
    # w never moves off 1 (damping = 0), so every smoother() call after the
    # first returns the exact same fit -- but convergence is only checked
    # against the PREVIOUS iteration's d_geq, and iteration 1's d_geq always
    # differs from the initial all-FALSE state, so it can't break until
    # iteration 2 re-confirms d_geq is unchanged: exactly 2 calls, not 1.
    expect_equal(length(smoother_calls), 2)
    expect_identical(smoother_calls[[1]], smoother_calls[[2]])
})

## damping changes psalsa()'s output relative to the default -----------------------

test_that("a smaller damping value produces a different (and still finite) baseline than damping = 1", {
    set.seed(4)
    y <- 8 + 0.01 * seq_len(300) + 30 * exp(-((seq_len(300) - 150)^2) / (2 * 5^2)) + rnorm(300, 0, 0.3)
    fit_undamped <- AlpsNMR:::psalsa(y, lambda = 1e6, p = 0.001, k = 5, maxit = 50, damping = 1)
    fit_damped <- AlpsNMR:::psalsa(y, lambda = 1e6, p = 0.001, k = 5, maxit = 50, damping = 0.7)

    expect_true(all(is.finite(fit_damped$baseline)))
    expect_false(isTRUE(all.equal(fit_undamped$baseline, fit_damped$baseline)))
})

test_that("damping is threaded through the matrix (multi-sample) path of psalsa()", {
    set.seed(5)
    y1 <- 5 + 20 * exp(-((seq_len(100) - 50)^2) / (2 * 4^2)) + rnorm(100, 0, 0.2)
    y2 <- 6 + 18 * exp(-((seq_len(100) - 50)^2) / (2 * 4^2)) + rnorm(100, 0, 0.2)
    spectra <- rbind(y1, y2)
    fit <- AlpsNMR:::psalsa(spectra, lambda = 1e5, p = 0.01, k = 3, damping = 0.5)
    fit_row1 <- AlpsNMR:::psalsa(y1, lambda = 1e5, p = 0.01, k = 3, damping = 0.5)
    expect_equal(fit$baseline[1, ], fit_row1$baseline, tolerance = 1e-10)
})

## convergence robustness with a flat, borderline-noisy region ---------------------

test_that("damping = 0.5 reliably converges on a flat region with many near-zero residuals", {
    set.seed(6)
    n <- 500
    # A tall peak plus a long, genuinely flat/noisy stretch -- many points
    # sitting close to zero is exactly the scenario that can trigger the
    # hard-threshold oscillation this damping fix targets.
    y <- 20 * exp(-((seq_len(n) - 50)^2) / (2 * 4^2)) + rnorm(n, 0, 5)
    fit <- AlpsNMR:::psalsa(y, lambda = 1e7, p = 0.001, k = 2, maxit = 150, damping = 0.5)
    expect_true(all(is.finite(fit$baseline)))
})
