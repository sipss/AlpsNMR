## diff2_penalty_weighted ------------------------------------------------------

test_that("diff2_penalty_weighted with a constant lambda matches lambda * diff2_penalty", {
    n <- 50
    lambda <- 3.7
    old <- lambda * AlpsNMR:::diff2_penalty(n)
    new_scalar <- AlpsNMR:::diff2_penalty_weighted(n, lambda)
    new_full_profile <- AlpsNMR:::diff2_penalty_weighted(n, rep(lambda, n))
    expect_equal(as.matrix(old), as.matrix(new_scalar), tolerance = 1e-10)
    expect_equal(as.matrix(old), as.matrix(new_full_profile), tolerance = 1e-10)
})

test_that("diff2_penalty_weighted with a position-varying lambda differs from any constant-lambda penalty", {
    n <- 30
    lambda_vec <- c(rep(1, 15), rep(100, 15))
    weighted <- AlpsNMR:::diff2_penalty_weighted(n, lambda_vec)
    constant <- AlpsNMR:::diff2_penalty_weighted(n, 1)
    expect_false(isTRUE(all.equal(as.matrix(weighted), as.matrix(constant))))
    # still symmetric, still a valid penalty:
    expect_equal(as.matrix(weighted), t(as.matrix(weighted)), tolerance = 1e-10)
})

## psalsa_core: vector p/k -----------------------------------------------------

test_that("psalsa_core with vector p/k matches scalar p/k when the vector is constant", {
    set.seed(1)
    y <- 5 + 20 * exp(-((seq_len(100) - 50)^2) / (2 * 4^2)) + rnorm(100, 0, 0.2)
    smoother <- function(y, w) AlpsNMR:::whit1d(y, 1e6 * AlpsNMR:::diff2_penalty(length(y)), w)

    scalar_fit <- AlpsNMR:::psalsa_core(y, smoother, p = 0.01, k = 3)
    vector_fit <- AlpsNMR:::psalsa_core(y, smoother, p = rep(0.01, 100), k = rep(3, 100))
    expect_equal(scalar_fit, vector_fit, tolerance = 1e-10)
})

test_that("psalsa_core's p[d_geq]/k[d_geq] indexing stays aligned (regression test)", {
    # A naive `p * exp(-d[d_geq] / k)` with vector p/k (no broadcasting to
    # full length first) would misalign p/k against the SUBSET d[d_geq]
    # rather than the point it actually belongs to. Use very different p/k
    # in the two halves and confirm the peak-side reweighting in each half
    # used ITS OWN half's k (a large k barely discounts a tall residual; a
    # tiny k discounts it almost to zero).
    set.seed(2)
    n <- 200
    y <- rep(1, n)
    y[50] <- 1 + 50 # a tall spike in the first half
    y[150] <- 1 + 50 # an equally tall spike in the second half
    k_vec <- c(rep(1e6, 100), rep(1e-6, 100)) # first half: barely discounts; second: nearly zeroes out
    smoother <- function(y, w) AlpsNMR:::whit1d(y, 1e2 * AlpsNMR:::diff2_penalty(n), w)

    # Run psalsa_core for a single iteration's worth of reweighting logic directly:
    p_vec <- rep(0.01, n)
    w <- rep(1, n)
    s <- smoother(y, w)
    d <- y - s
    d_geq <- d >= 0
    k_full <- rep_len(k_vec, n)
    p_full <- rep_len(p_vec, n)
    w[d_geq] <- p_full[d_geq] * exp(-d[d_geq] / k_full[d_geq])

    # The point at index 50 (huge k -> exp(~0) ~ 1 -> weight ~ p) should end
    # up with a much larger weight than index 150 (tiny k -> exp(-huge) ~ 0):
    expect_true(w[50] > w[150] * 100)
})

## psalsa_core: k defaulting still scalar-only ----------------------------------

test_that("psalsa_core's k = -1 auto-default only triggers for a single value", {
    set.seed(3)
    y <- 5 + 10 * exp(-((seq_len(80) - 40)^2) / (2 * 3^2)) + rnorm(80, 0, 0.1)
    smoother <- function(y, w) AlpsNMR:::whit1d(y, 1e5 * AlpsNMR:::diff2_penalty(length(y)), w)
    # A vector k containing -1 should NOT be auto-defaulted (used literally, floored):
    result <- AlpsNMR:::psalsa_core(y, smoother, p = 0.01, k = rep(-1, 80))
    expect_true(all(is.finite(result)))
})

## psalsa(): end-to-end with position-varying lambda/p/k -----------------------

test_that("psalsa() with vector lambda/p/k produces a finite baseline the right length", {
    set.seed(4)
    n <- 300
    y <- 8 + 0.01 * seq_len(n) + 30 * exp(-((seq_len(n) - 150)^2) / (2 * 5^2)) + rnorm(n, 0, 0.3)
    lambda_vec <- seq(1e4, 1e8, length.out = n)
    p_vec <- rep(0.001, n)
    k_vec <- rep(5, n)

    fit <- AlpsNMR:::psalsa(y, lambda = lambda_vec, p = p_vec, k = k_vec)
    expect_equal(length(fit$baseline), n)
    expect_true(all(is.finite(fit$baseline)))
    expect_equal(fit$corrected, y - fit$baseline)
})

test_that("psalsa() applies the same position profile to every row of a matrix", {
    set.seed(5)
    n <- 100
    y1 <- 5 + 10 * exp(-((seq_len(n) - 50)^2) / (2 * 3^2)) + rnorm(n, 0, 0.1)
    y2 <- 7 + 8 * exp(-((seq_len(n) - 50)^2) / (2 * 3^2)) + rnorm(n, 0, 0.1)
    spectra <- rbind(y1, y2)
    lambda_vec <- c(rep(1e3, 50), rep(1e7, 50))

    fit <- AlpsNMR:::psalsa(spectra, lambda = lambda_vec, p = 0.001, k = 3)
    expect_equal(dim(fit$baseline), dim(spectra))
    fit1 <- AlpsNMR:::psalsa(y1, lambda = lambda_vec, p = 0.001, k = 3)
    expect_equal(fit$baseline[1, ], fit1$baseline, tolerance = 1e-10)
})

## spatial_profile_1d -----------------------------------------------------------

test_that("spatial_profile_1d's log transform stays positive across knots of wildly different scale", {
    profile <- AlpsNMR:::spatial_profile_1d(
        n = 100, frac_mid = c(0.1, 0.5, 0.9), values = c(1e3, 1e9, 1e4), transform = "log"
    )
    expect_equal(length(profile), 100)
    expect_true(all(profile > 0))
    expect_true(all(is.finite(profile)))
})

test_that("spatial_profile_1d's logit transform stays within (0, p_max)", {
    profile <- AlpsNMR:::spatial_profile_1d(
        n = 100, frac_mid = c(0.1, 0.5, 0.9), values = c(0.001, 0.04, 0.002),
        transform = "logit", p_max = 0.05
    )
    expect_true(all(profile > 0 & profile < 0.05))
})

test_that("spatial_profile_1d falls back to a constant profile with fewer than 2 knots", {
    profile <- AlpsNMR:::spatial_profile_1d(n = 50, frac_mid = 0.5, values = 1e6, transform = "log")
    expect_equal(profile, rep(1e6, 50), tolerance = 1e-8)
})

test_that("spatial_profile_1d interpolates smoothly near its knots (roughly recovers knot values)", {
    profile <- AlpsNMR:::spatial_profile_1d(
        n = 1000, frac_mid = c(0.2, 0.8), values = c(1e5, 1e8), transform = "log"
    )
    expect_equal(profile[200], 1e5, tolerance = 0.05)
    expect_equal(profile[800], 1e8, tolerance = 0.05)
    expect_true(profile[500] > profile[200] && profile[500] < profile[800])
})

## resample_profile_1d -----------------------------------------------------------

test_that("resample_profile_1d is a no-op when lengths already match", {
    profile <- seq_len(50)
    expect_identical(AlpsNMR:::resample_profile_1d(profile, 50), profile)
})

test_that("resample_profile_1d preserves overall shape when stretching/shrinking", {
    profile <- c(rep(1, 50), rep(10, 50)) # low then high, split at the midpoint
    stretched <- AlpsNMR:::resample_profile_1d(profile, 200)
    expect_equal(length(stretched), 200)
    expect_true(mean(stretched[1:90]) < mean(stretched[111:200]))

    shrunk <- AlpsNMR:::resample_profile_1d(profile, 20)
    expect_equal(length(shrunk), 20)
    expect_true(mean(shrunk[1:9]) < mean(shrunk[12:20]))
})

## tune_psalsa_region_params_1d --------------------------------------------------

test_that("tune_psalsa_region_params_1d skips regions with no observed peaks", {
    pooled_regions <- list(
        csnr = 0.03,
        regions = data.frame(
            region = 1:3, frac_lo = c(0, 1 / 3, 2 / 3), frac_hi = c(1 / 3, 2 / 3, 1),
            density = c(0.02, 0, 0.015), fwhm_q1 = c(5, NA, 6), fwhm_q3 = c(10, NA, 12)
        )
    )
    # min_peaks = 0 disables adjacent-region merging (see
    # test-psalsa_region_anchor.R for merging's own tests), isolating the
    # "skip an empty region" behaviour this test targets.
    result <- AlpsNMR:::tune_psalsa_region_params_1d(
        pooled_regions, n_total = 900, peak_shape = "gaussian", n_synthetic = 3, optim_maxit = 20,
        min_peaks = 0
    )
    expect_equal(sort(result$region), c(1, 3))
})

test_that("tune_psalsa_region_params_1d returns NULL when no region has observed peaks", {
    pooled_regions <- list(
        csnr = 0.03,
        regions = data.frame(
            region = 1:2, frac_lo = c(0, 0.5), frac_hi = c(0.5, 1),
            density = c(0, 0), fwhm_q1 = c(NA, NA), fwhm_q3 = c(NA, NA)
        )
    )
    result <- AlpsNMR:::tune_psalsa_region_params_1d(
        pooled_regions, n_total = 200, peak_shape = "gaussian", n_synthetic = 3, optim_maxit = 20
    )
    expect_null(result)
})

## tune_psalsa_spatial(): end-to-end ----------------------------------------------

test_that("tune_psalsa_spatial runs end-to-end and returns position-varying profiles", {
    set.seed(6)
    n <- 500
    baseline <- 10 + 3 * sin(seq_len(n) / 100)
    # Crowded region (many small peaks) in the first half, sparse in the second:
    crowded <- Reduce(`+`, lapply(seq(40, 200, by = 20), function(ctr) {
        6 * exp(-((seq_len(n) - ctr)^2) / (2 * 3^2))
    }))
    sparse <- 40 * exp(-((seq_len(n) - 400)^2) / (2 * 6^2))
    y <- baseline + crowded + sparse + rnorm(n, 0, 0.4)

    # min_peaks = 0 disables adjacent-region merging so this fixture's own
    # 8-region split is what gets tuned -- merging's effect on the final
    # profile is covered separately in test-psalsa_region_anchor.R.
    result <- AlpsNMR:::tune_psalsa_spatial(y, num_regions = 8, n_synthetic = 3, optim_maxit = 30, min_peaks = 0)
    expect_equal(length(result$baseline), n)
    expect_true(all(is.finite(result$baseline)))
    expect_equal(length(result$lambda), n)
    expect_equal(length(result$p), n)
    expect_equal(length(result$k), n)
    expect_true(all(result$lambda > 0))
    expect_true(all(result$p > 0 & result$p < 0.05))
    expect_true(all(result$k > 0))
    expect_true(is.data.frame(result$region_params))
    expect_true(nrow(result$region_params) >= 2)
})

test_that("tune_psalsa_spatial falls back to tune_psalsa when fewer than 2 regions have signal", {
    set.seed(7)
    n <- 200
    # A single narrow burst of peaks confined to a tiny sliver -- with many
    # regions requested, at most 1 will ever see a peak.
    y <- rep(5, n) + rnorm(n, 0, 0.05)
    y[100] <- y[100] + 30
    y[101] <- y[101] + 28

    result <- AlpsNMR:::tune_psalsa_spatial(y, num_regions = 50, n_synthetic = 3, optim_maxit = 20)
    # Falls back to tune_psalsa()'s output shape: scalar lambda/p/k, no region_params.
    expect_equal(length(result$baseline), n)
    expect_true(is.null(result$region_params))
    expect_true(length(result$lambda) == 1)
})

test_that("tune_psalsa_spatial handles a list of differently-lengthed spectra", {
    set.seed(8)
    make_y <- function(n) {
        baseline <- 8 + 2 * sin(seq_len(n) / 50)
        peaks <- Reduce(`+`, lapply(seq(0.1, 0.9, by = 0.2) * n, function(ctr) {
            10 * exp(-((seq_len(n) - ctr)^2) / (2 * 4^2))
        }))
        baseline + peaks + rnorm(n, 0, 0.3)
    }
    y_list <- list(make_y(400), make_y(420))

    result <- AlpsNMR:::tune_psalsa_spatial(y_list, num_regions = 6, n_synthetic = 3, optim_maxit = 20)
    expect_equal(length(result$baseline), 2)
    expect_equal(length(result$baseline[[1]]), 400)
    expect_equal(length(result$baseline[[2]]), 420)
    expect_true(all(is.finite(result$baseline[[1]])))
    expect_true(all(is.finite(result$baseline[[2]])))
})
