## analyze_signal_regions_1d: per-region amplitude ------------------------------

test_that("analyze_signal_regions_1d reports a much larger amplitude for a region with taller peaks", {
    set.seed(1)
    n <- 2000
    y <- rep(1, n) + rnorm(n, 0, 0.05)
    # Tall peaks in the first quarter, small peaks in the third quarter:
    for (ctr in seq(50, 450, by = 60)) y <- y + 5000 * exp(-((seq_len(n) - ctr)^2) / (2 * 3^2))
    for (ctr in seq(1050, 1450, by = 60)) y <- y + 50 * exp(-((seq_len(n) - ctr)^2) / (2 * 3^2))

    result <- analyze_signal_regions_1d(y, num_regions = 4)
    tall_region <- result$regions[1, ]
    small_region <- result$regions[3, ]
    expect_true(!is.na(tall_region$amplitude_q2) && !is.na(small_region$amplitude_q2))
    expect_true(tall_region$amplitude_q2 > 10 * small_region$amplitude_q2)
})

test_that("analyze_signal_regions_1d gives NA amplitude for a region with no detected peaks", {
    set.seed(2)
    n <- 1000
    y <- rep(1, n) + rnorm(n, 0, 0.05) # no peaks anywhere
    result <- analyze_signal_regions_1d(y, num_regions = 4)
    expect_true(all(is.na(result$regions$amplitude_q1)))
    expect_true(all(is.na(result$regions$amplitude_q2)))
    expect_true(all(is.na(result$regions$amplitude_q3)))
    expect_true(all(result$regions$density == 0))
})

## pool_signal_stats_regions_1d: pooling amplitude -------------------------------

test_that("pool_signal_stats_regions_1d pools amplitude (median) across signals", {
    set.seed(3)
    make_signal <- function(n, height) {
        y <- rep(1, n) + rnorm(n, 0, 0.05)
        for (ctr in seq(50, 450, by = 60)) y <- y + height * exp(-((seq_len(n) - ctr)^2) / (2 * 3^2))
        y
    }
    y_list <- list(make_signal(1000, 100), make_signal(1000, 300), make_signal(1000, 200))
    pooled <- pool_signal_stats_regions_1d(y_list, num_regions = 4)
    expect_true(all(c("amplitude_q1", "amplitude_q2", "amplitude_q3") %in% colnames(pooled$regions)))
    expect_true(!is.na(pooled$regions$amplitude_q2[1]))
    expect_true(pooled$regions$amplitude_q2[1] > 50) # roughly tracks 0.35^-1 * median(100,200,300)-ish scale
    expect_true(pooled$regions$amplitude_q1[1] <= pooled$regions$amplitude_q2[1])
    expect_true(pooled$regions$amplitude_q2[1] <= pooled$regions$amplitude_q3[1])
})

## merge_sparse_regions_1d: amplitude recombination ------------------------------

test_that("merge_sparse_regions_1d weights amplitude by each member's own peak count", {
    regions <- data.frame(
        region = 1:2, frac_lo = c(0, 0.5), frac_hi = c(0.5, 1),
        density_q1 = c(0.001, 0.02), density_q2 = c(0.001, 0.02), density_q3 = c(0.001, 0.02),
        fwhm_q1 = c(5, 5), fwhm_q3 = c(10, 10),
        amplitude_q1 = c(1000, 100), amplitude_q2 = c(1000, 100), amplitude_q3 = c(1000, 100)
    )
    # region 2 has ~20x the estimated peak count of region 1 (given equal
    # span), so the merged amplitude should sit much closer to 100 than 1000:
    merged <- merge_sparse_regions_1d(regions, n_total = 1000, min_peaks = 1)
    expect_equal(nrow(merged), 1)
    expect_true(merged$amplitude_q2 < 200)
    expect_true(merged$amplitude_q2 > 100)
})

test_that("merge_sparse_regions_1d defaults amplitude to NA when the input regions table lacks it", {
    regions <- data.frame(
        region = 1:2, frac_lo = c(0, 0.5), frac_hi = c(0.5, 1),
        density_q1 = c(0.1, 0.1), density_q2 = c(0.1, 0.1), density_q3 = c(0.1, 0.1),
        fwhm_q1 = c(5, 5), fwhm_q3 = c(10, 10)
    )
    merged <- merge_sparse_regions_1d(regions, n_total = 1000, min_peaks = 1)
    expect_true(all(c("amplitude_q1", "amplitude_q2", "amplitude_q3") %in% colnames(merged)))
    expect_true(all(is.na(merged$amplitude_q1)))
    expect_true(all(is.na(merged$amplitude_q2)))
    expect_true(all(is.na(merged$amplitude_q3)))
})

## gen_synthetic_1d_regions: region-specific amplitude ---------------------------

test_that("gen_synthetic_1d_regions uses region_profile's amplitude_q1/q2/q3 when present", {
    region_profile <- data.frame(
        region = 1:2, frac_lo = c(0, 0.5), frac_hi = c(0.5, 1),
        density_q1 = c(0.05, 0.05), density_q2 = c(0.05, 0.05), density_q3 = c(0.05, 0.05),
        fwhm_q1 = c(5, 5), fwhm_q3 = c(10, 10),
        amplitude_q1 = c(10000, 10), amplitude_q2 = c(10000, 10), amplitude_q3 = c(10000, 10)
    )
    sig <- gen_synthetic_1d_regions(n = 1000, region_profile = region_profile, csnr = 0, seed = 1)
    tall_heights <- sig$peak_info$height[sig$peak_info$center <= 500]
    small_heights <- sig$peak_info$height[sig$peak_info$center > 500]
    expect_true(length(tall_heights) > 0 && length(small_heights) > 0)
    expect_true(median(tall_heights) > 100 * median(small_heights))
})

test_that("gen_synthetic_1d_regions falls back to A when amplitude is absent (unchanged behaviour)", {
    region_profile <- data.frame(
        region = 1, frac_lo = 0, frac_hi = 1,
        density_q1 = 0.02, density_q2 = 0.02, density_q3 = 0.02, fwhm_q1 = 5, fwhm_q3 = 10
    )
    set.seed(5)
    sig_no_amp <- gen_synthetic_1d_regions(n = 500, region_profile = region_profile, A = 7, seed = 1)
    set.seed(5)
    region_profile_na <- cbind(region_profile, amplitude_q1 = NA_real_, amplitude_q2 = NA_real_, amplitude_q3 = NA_real_)
    sig_na_amp <- gen_synthetic_1d_regions(n = 500, region_profile = region_profile_na, A = 7, seed = 1)
    expect_equal(sig_no_amp$peak_info$height, sig_na_amp$peak_info$height)
})

## tune_psalsa_region_params_1d: amplitude fallback layers -----------------------

test_that("tune_psalsa_region_params_1d falls back to the ORIGINAL behaviour when neither amplitude nor noise_sd is available", {
    pooled_regions <- list(
        csnr = 0.03,
        regions = data.frame(
            region = 1, frac_lo = 0, frac_hi = 1,
            density_q1 = 0.02, density_q2 = 0.02, density_q3 = 0.02, fwhm_q1 = 5, fwhm_q3 = 10
        )
    )
    set.seed(6)
    result <- tune_psalsa_region_params_1d(
        pooled_regions, n_total = 400, peak_shape = "gaussian", n_synthetic = 3, optim_maxit = 20,
        min_peaks = 0
    )
    expect_equal(nrow(result), 1)
    expect_true(is.finite(result$lambda))
})

test_that("tune_psalsa_region_params_1d uses signal-wide amplitude (noise_sd/csnr) when a region's own amplitude is NA", {
    pooled_regions <- list(
        csnr = 0.03, noise_sd = 3,
        regions = data.frame(
            region = 1, frac_lo = 0, frac_hi = 1,
            density_q1 = 0.02, density_q2 = 0.02, density_q3 = 0.02, fwhm_q1 = 5, fwhm_q3 = 10
        )
    )
    set.seed(7)
    result <- tune_psalsa_region_params_1d(
        pooled_regions, n_total = 400, peak_shape = "gaussian", n_synthetic = 3, optim_maxit = 20,
        min_peaks = 0
    )
    expect_equal(nrow(result), 1)
    expect_true(is.finite(result$lambda))
})

test_that("a region's amplitude materially changes its tuned parameters (not just plumbed through inertly)", {
    # Same density/fwhm/seed, only amplitude differs. The exact direction
    # k_mult moves isn't a simple monotonic function of amplitude alone (it
    # emerges from the whole area-recovery objective against each region's
    # own resulting csnr = noise_sd/amplitude) -- what this actually checks
    # is that amplitude is reaching the tuning search at all and changing
    # its outcome, which a purely-plumbing bug (e.g. amplitude silently
    # ignored) would not produce.
    pooled_regions_tall <- list(
        noise_sd = 1, csnr = 0.001,
        regions = data.frame(
            region = 1, frac_lo = 0, frac_hi = 1,
            density_q1 = 0.02, density_q2 = 0.02, density_q3 = 0.02, fwhm_q1 = 5, fwhm_q3 = 10,
            amplitude_q1 = 100000, amplitude_q2 = 100000, amplitude_q3 = 100000
        )
    )
    pooled_regions_small <- list(
        noise_sd = 1, csnr = 0.001,
        regions = data.frame(
            region = 1, frac_lo = 0, frac_hi = 1,
            density_q1 = 0.02, density_q2 = 0.02, density_q3 = 0.02, fwhm_q1 = 5, fwhm_q3 = 10,
            amplitude_q1 = 10, amplitude_q2 = 10, amplitude_q3 = 10
        )
    )
    set.seed(8)
    result_tall <- tune_psalsa_region_params_1d(
        pooled_regions_tall, n_total = 500, peak_shape = "gaussian", n_synthetic = 3, optim_maxit = 30,
        min_peaks = 0
    )
    set.seed(8)
    result_small <- tune_psalsa_region_params_1d(
        pooled_regions_small, n_total = 500, peak_shape = "gaussian", n_synthetic = 3, optim_maxit = 30,
        min_peaks = 0
    )
    expect_false(isTRUE(all.equal(result_tall$lambda, result_small$lambda)))
    expect_false(isTRUE(all.equal(result_tall$k_mult, result_small$k_mult)))
})
