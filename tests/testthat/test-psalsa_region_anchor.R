## tune_psalsa_region_params_1d: theta0/p_max anchoring -------------------------

test_that("tune_psalsa_region_params_1d with theta0 = NULL keeps the previous, unanchored behaviour", {
    pooled_regions <- list(
        csnr = 0.03,
        regions = data.frame(
            region = 1:2, frac_lo = c(0, 0.5), frac_hi = c(0.5, 1),
            density = c(0.02, 0.02), fwhm_q1 = c(5, 5), fwhm_q3 = c(10, 10)
        )
    )
    set.seed(1)
    # min_peaks = 0 disables adjacent-region merging (its own tests are
    # below), keeping this fixture's 2 regions separate.
    result <- AlpsNMR:::tune_psalsa_region_params_1d(
        pooled_regions, n_total = 400, peak_shape = "gaussian", n_synthetic = 3, optim_maxit = 30,
        min_peaks = 0
    )
    expect_equal(nrow(result), 2)
    expect_true(all(is.finite(result$lambda)))
})

test_that("a custom theta0 anchor pulls a sparse, weakly-constrained region's tuned value toward it", {
    # A single region with almost no peaks: its own objective is close to
    # flat, so the regularization term should dominate and pin the result
    # near whichever anchor is supplied.
    pooled_regions <- list(
        csnr = 0.03,
        regions = data.frame(
            region = 1, frac_lo = 0, frac_hi = 1, density = 0.0005, fwhm_q1 = 5, fwhm_q3 = 10
        )
    )
    p_max <- 0.05
    anchor_lo <- c(log(1e4), stats::qlogis(0.001 / p_max), log(5))
    anchor_hi <- c(log(1e9), stats::qlogis(0.001 / p_max), log(5))

    set.seed(10)
    result_lo <- AlpsNMR:::tune_psalsa_region_params_1d(
        pooled_regions, n_total = 400, peak_shape = "gaussian", n_synthetic = 3, optim_maxit = 60,
        theta0 = anchor_lo, p_max = p_max
    )
    set.seed(10)
    result_hi <- AlpsNMR:::tune_psalsa_region_params_1d(
        pooled_regions, n_total = 400, peak_shape = "gaussian", n_synthetic = 3, optim_maxit = 60,
        theta0 = anchor_hi, p_max = p_max
    )

    # Same sparse region, same seed, only the anchor differs -- the tuned
    # lambda should land closer to ITS OWN anchor than to the other one:
    expect_true(abs(log(result_lo$lambda) - log(1e4)) < abs(log(result_lo$lambda) - log(1e9)))
    expect_true(abs(log(result_hi$lambda) - log(1e9)) < abs(log(result_hi$lambda) - log(1e4)))
})

test_that("p_max is threaded through to the per-region call alongside theta0", {
    pooled_regions <- list(
        csnr = 0.03,
        regions = data.frame(
            region = 1, frac_lo = 0, frac_hi = 1, density = 0.02, fwhm_q1 = 5, fwhm_q3 = 10
        )
    )
    p_max_custom <- 0.2
    anchor <- c(log(1e6), stats::qlogis(0.05 / p_max_custom), log(10))
    set.seed(2)
    result <- AlpsNMR:::tune_psalsa_region_params_1d(
        pooled_regions, n_total = 300, peak_shape = "gaussian", n_synthetic = 3, optim_maxit = 30,
        theta0 = anchor, p_max = p_max_custom
    )
    # p stays a valid probability regardless of the non-default p_max:
    expect_true(result$p > 0 && result$p < p_max_custom)
})

test_that("region_lambda_k_weight/region_p_weight keep a very sparse region from collapsing (regression test)", {
    # A region with only 3 peaks in its own synthetic pool -- reproduces the
    # real failure mode found on the MTBLS242 dataset: at the GLOBAL search's
    # own (much weaker) regularization weights, this collapses to lambda
    # orders of magnitude below the anchor even though the anchor is passed.
    p_max <- 0.05
    anchor <- c(log(2.8e7), stats::qlogis(0.0005 / p_max), log(190))
    set.seed(55)
    y_sparse <- gen_synthetic_1d(n = 356, density = 0.003, fwhm_range = c(3, 3), csnr = 0.03, seed = 55)
    pool_r <- list(y_sparse)

    weak <- tune_psalsa_params_1d(pool_r, theta0 = anchor, p_max = p_max, lambda_k_weight = 0.02, p_weight = 0.4)
    strong <- tune_psalsa_params_1d(pool_r, theta0 = anchor, p_max = p_max, lambda_k_weight = 2, p_weight = 10)

    # The stronger, per-region defaults must land closer to the anchor (in
    # log-lambda space) than the weak, whole-spectrum-calibrated weights do:
    expect_true(abs(log(strong$lambda) - log(2.8e7)) < abs(log(weak$lambda) - log(2.8e7)))
})

test_that("tune_psalsa_spatial's default region weights keep every region's lambda within a bounded range of the pack", {
    set.seed(56)
    n <- 1500
    baseline <- 10 + 3 * sin(seq_len(n) / 200)
    # Dense peaks everywhere except one very sparse, narrow sliver:
    peaks <- Reduce(`+`, lapply(seq(20, 1450, by = 30), function(ctr) {
        if (ctr > 1300) return(rep(0, n))
        7 * exp(-((seq_len(n) - ctr)^2) / (2 * 3^2))
    }))
    y <- baseline + peaks + rnorm(n, 0, 0.3)

    result <- AlpsNMR:::tune_psalsa_spatial(y, num_regions = 15, n_synthetic = 3, optim_maxit = 40)
    log_lambda <- log(result$region_params$lambda)
    # No region's own lambda should be more than ~3 orders of magnitude from
    # the pack's median -- the collapse this regularization is meant to stop
    # showed up as a 4-5+ order-of-magnitude swing on real data:
    expect_true(max(abs(log_lambda - stats::median(log_lambda))) < log(1e3))
})

## merge_sparse_regions_1d ------------------------------------------------------

test_that("merge_sparse_regions_1d leaves regions alone when each already meets min_peaks", {
    regions <- data.frame(
        region = 1:3, frac_lo = c(0, 1 / 3, 2 / 3), frac_hi = c(1 / 3, 2 / 3, 1),
        density = c(0.1, 0.1, 0.1), fwhm_q1 = c(5, 5, 5), fwhm_q3 = c(10, 10, 10)
    )
    merged <- AlpsNMR:::merge_sparse_regions_1d(regions, n_total = 900, min_peaks = 15)
    # Each region alone has 0.1 * 300 = 30 estimated peaks, well above 15:
    expect_equal(nrow(merged), 3)
    expect_equal(merged$frac_lo, regions$frac_lo)
    expect_equal(merged$frac_hi, regions$frac_hi)
})

test_that("merge_sparse_regions_1d merges adjacent sparse regions until min_peaks is reached", {
    regions <- data.frame(
        region = 1:4, frac_lo = c(0, 0.25, 0.5, 0.75), frac_hi = c(0.25, 0.5, 0.75, 1),
        density = c(0.002, 0.002, 0.002, 0.002), fwhm_q1 = rep(5, 4), fwhm_q3 = rep(10, 4)
    )
    # Each region alone: 0.002 * 250 = 0.5 estimated peaks -- all 4 must
    # merge into one group to reach min_peaks = 2 (2 estimated peaks total,
    # only reached once every region has been folded in):
    merged <- AlpsNMR:::merge_sparse_regions_1d(regions, n_total = 1000, min_peaks = 2)
    expect_equal(nrow(merged), 1)
    expect_equal(merged$frac_lo, 0)
    expect_equal(merged$frac_hi, 1)
})

test_that("merge_sparse_regions_1d merges a trailing under-informed group backward", {
    regions <- data.frame(
        region = 1:3, frac_lo = c(0, 1 / 3, 2 / 3), frac_hi = c(1 / 3, 2 / 3, 1),
        density = c(0.1, 0.1, 0.0001), fwhm_q1 = c(5, 5, 5), fwhm_q3 = c(10, 10, 10)
    )
    # Region 1 alone reaches min_peaks (30 peaks); region 2 alone also
    # reaches it (30 peaks); region 3 alone falls far short (0.03 peaks) and
    # must merge backward into region 2 rather than stay its own group:
    merged <- AlpsNMR:::merge_sparse_regions_1d(regions, n_total = 900, min_peaks = 15)
    expect_equal(nrow(merged), 2)
    expect_equal(merged$frac_lo, c(0, 1 / 3))
    expect_equal(merged$frac_hi, c(1 / 3, 1))
})

test_that("merge_sparse_regions_1d doesn't force-merge a trailing empty region into an already well-informed neighbour", {
    regions <- data.frame(
        region = 1:2, frac_lo = c(0, 0.5), frac_hi = c(0.5, 1),
        density = c(0.2, 0), fwhm_q1 = c(5, NA), fwhm_q3 = c(10, NA)
    )
    merged <- AlpsNMR:::merge_sparse_regions_1d(regions, n_total = 200, min_peaks = 15)
    # Region 1 alone reaches min_peaks (20 peaks) and closes its own group;
    # region 2 is genuinely empty (0 estimated peaks) and gains nothing from
    # merging, so it's left as its own (skippable, density-0) group rather
    # than diluting region 1's density and span:
    expect_equal(nrow(merged), 2)
    expect_equal(merged$density[2], 0)
    expect_true(is.na(merged$fwhm_q1[2]))
})

test_that("merge_sparse_regions_1d keeps an all-empty spectrum as a single density-0 group", {
    regions <- data.frame(
        region = 1:3, frac_lo = c(0, 1 / 3, 2 / 3), frac_hi = c(1 / 3, 2 / 3, 1),
        density = c(0, 0, 0), fwhm_q1 = c(NA, NA, NA), fwhm_q3 = c(NA, NA, NA)
    )
    merged <- AlpsNMR:::merge_sparse_regions_1d(regions, n_total = 900, min_peaks = 15)
    expect_equal(nrow(merged), 1)
    expect_equal(merged$density, 0)
    expect_true(is.na(merged$fwhm_q1))
})

## tune_psalsa_region_params_1d: min_peaks integration ---------------------------

test_that("min_peaks controls whether sparse adjacent regions get merged before tuning", {
    pooled_regions <- list(
        csnr = 0.03,
        regions = data.frame(
            region = 1:2, frac_lo = c(0, 0.5), frac_hi = c(0.5, 1),
            density = c(0.001, 0.001), fwhm_q1 = c(5, 5), fwhm_q3 = c(10, 10)
        )
    )
    set.seed(20)
    unmerged <- AlpsNMR:::tune_psalsa_region_params_1d(
        pooled_regions, n_total = 400, peak_shape = "gaussian", n_synthetic = 3, optim_maxit = 20,
        min_peaks = 0
    )
    set.seed(20)
    merged <- AlpsNMR:::tune_psalsa_region_params_1d(
        pooled_regions, n_total = 400, peak_shape = "gaussian", n_synthetic = 3, optim_maxit = 20,
        min_peaks = 15
    )
    expect_equal(nrow(unmerged), 2)
    expect_equal(nrow(merged), 1)
})

test_that("with adjacent-region merging (default min_peaks), a sparse region's tuned value stays close to the anchor without needing the strongest regularization weights", {
    # The same sparse-region scenario that originally motivated
    # region_lambda_k_weight/region_p_weight: with min_peaks doing its job,
    # even a MUCH weaker regularization weight (matching the whole-spectrum
    # defaults) should no longer see the sparse region collapse, since by
    # the time it's tuned it's been merged into a group with real data.
    pooled_regions <- list(
        csnr = 0.03,
        regions = data.frame(
            region = 1:2, frac_lo = c(0, 0.5), frac_hi = c(0.5, 1),
            density = c(0.05, 0.0005), fwhm_q1 = c(5, 5), fwhm_q3 = c(10, 10)
        )
    )
    p_max <- 0.05
    anchor <- c(log(2.8e7), stats::qlogis(0.0005 / p_max), log(190))
    set.seed(30)
    result <- AlpsNMR:::tune_psalsa_region_params_1d(
        pooled_regions, n_total = 2000, peak_shape = "gaussian", n_synthetic = 3, optim_maxit = 40,
        theta0 = anchor, p_max = p_max,
        region_lambda_k_weight = 0.02, region_p_weight = 0.4, # whole-spectrum-strength weights
        min_peaks = 15
    )
    # Both original regions merge into one well-populated group -- no
    # collapse, and no need for the stronger weights to prevent it:
    expect_equal(nrow(result), 1)
    expect_true(abs(log(result$lambda) - log(2.8e7)) < log(100))
})

## tune_psalsa_spatial(): anchor wiring end-to-end -------------------------------

test_that("tune_psalsa_spatial's per-region tuning stays close to the spectrum-wide anchor for a sparse region", {
    set.seed(11)
    n <- 1000
    baseline <- 10 + 3 * sin(seq_len(n) / 150)
    # Plenty of peaks everywhere EXCEPT one narrow, deliberately sparse stretch:
    peaks <- Reduce(`+`, lapply(seq(20, 950, by = 40), function(ctr) {
        if (ctr > 400 && ctr < 600) {
            return(rep(0, n))
        }
        8 * exp(-((seq_len(n) - ctr)^2) / (2 * 3^2))
    }))
    y <- baseline + peaks + rnorm(n, 0, 0.3)

    # min_peaks = 0 disables adjacent-region merging here so the sparse
    # region this test targets stays its own knot instead of being merged
    # away -- merging's own effect is covered separately below.
    result <- AlpsNMR:::tune_psalsa_spatial(y, num_regions = 10, n_synthetic = 3, optim_maxit = 30, min_peaks = 0)
    expect_true(is.data.frame(result$region_params))
    # No single region's lambda should be many orders of magnitude away from
    # the pack's median -- the whole point of the anchor:
    log_lambda <- log(result$region_params$lambda)
    expect_true(max(abs(log_lambda - stats::median(log_lambda))) < log(1e4))
})
