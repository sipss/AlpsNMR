## compute_region_bounds_1d --------------------------------------------------

test_that("compute_region_bounds_1d partitions 1:n contiguously with no gaps or overlaps", {
    bounds <- compute_region_bounds_1d(100, 20)
    expect_equal(nrow(bounds), 20)
    expect_equal(bounds$lo[1], 1)
    expect_equal(bounds$hi[nrow(bounds)], 100)
    expect_equal(bounds$lo[-1], bounds$hi[-nrow(bounds)] + 1)
    expect_true(all(bounds$hi >= bounds$lo))
})

test_that("compute_region_bounds_1d caps num_regions at n", {
    bounds <- compute_region_bounds_1d(5, 20)
    expect_equal(nrow(bounds), 5)
    expect_equal(bounds$lo, 1:5)
    expect_equal(bounds$hi, 1:5)
})

## analyze_signal_regions_1d --------------------------------------------------

test_that("analyze_signal_regions_1d reports higher density in the region actually containing peaks", {
    set.seed(1)
    n <- 2000
    y <- rep(1, n)
    # A cluster of tall, narrow peaks confined to the first tenth of the signal:
    peak_centers <- seq(50, 150, by = 20)
    for (ctr in peak_centers) {
        y <- y + 20 * exp(-((seq_len(n) - ctr)^2) / (2 * 3^2))
    }
    y <- y + rnorm(n, 0, 0.05)

    result <- analyze_signal_regions_1d(y, num_regions = 20)
    expect_equal(nrow(result$regions), 20)
    expect_true(all(c("frac_lo", "frac_hi", "density", "fwhm_q1", "fwhm_q3") %in% colnames(result$regions)))

    # The peaks sit in roughly the first region (indices 50-150 of 2000 ~ region 1-2):
    busy_region_density <- result$regions$density[1]
    empty_region_density <- result$regions$density[15]
    expect_true(busy_region_density > empty_region_density)
    expect_equal(empty_region_density, 0)
})

## pool_signal_stats_regions_1d ------------------------------------------------

test_that("pool_signal_stats_regions_1d pools per-region stats across signals of different lengths", {
    set.seed(2)
    make_signal <- function(n, ctr_frac) {
        y <- rep(1, n) + rnorm(n, 0, 0.05)
        ctr <- round(ctr_frac * n)
        y <- y + 15 * exp(-((seq_len(n) - ctr)^2) / (2 * 3^2))
        y
    }
    # Same RELATIVE peak location (10% into the signal) at two different lengths:
    y_list <- list(make_signal(1000, 0.1), make_signal(1500, 0.1))

    pooled <- pool_signal_stats_regions_1d(y_list, num_regions = 10)
    expect_equal(nrow(pooled$regions), 10)
    expect_true(is.finite(pooled$noise_sd))
    expect_true(is.finite(pooled$csnr))
    # Region 1 (the first 10%) should show nonzero density; a far region should not:
    expect_true(pooled$regions$density_q2[1] > 0)
    expect_equal(pooled$regions$density_q2[8], 0)
})

## gen_synthetic_1d_regions -----------------------------------------------------

test_that("gen_synthetic_1d_regions places peaks only in regions with positive density", {
    region_profile <- data.frame(
        region = 1:2, frac_lo = c(0, 0.5), frac_hi = c(0.5, 1),
        density_q1 = c(0.05, 0), density_q2 = c(0.05, 0), density_q3 = c(0.05, 0),
        fwhm_q1 = c(5, NA), fwhm_q3 = c(10, NA)
    )
    sig <- gen_synthetic_1d_regions(n = 1000, region_profile = region_profile, seed = 1)

    expect_equal(length(sig$y), 1000)
    expect_true(nrow(sig$peak_info) > 0)
    # No peak center should land in the empty second half:
    expect_true(all(sig$peak_info$center <= 500))
})

test_that("gen_synthetic_1d_regions returns an empty peak_info when every region is empty", {
    region_profile <- data.frame(
        region = 1, frac_lo = 0, frac_hi = 1,
        density_q1 = 0, density_q2 = 0, density_q3 = 0, fwhm_q1 = NA_real_, fwhm_q3 = NA_real_
    )
    sig <- gen_synthetic_1d_regions(n = 200, region_profile = region_profile, seed = 1)
    expect_equal(nrow(sig$peak_info), 0)
    expect_equal(length(sig$y), 200)
})

test_that("gen_synthetic_1d_regions reproduces denser regions with more peaks than sparser ones", {
    region_profile <- data.frame(
        region = 1:2, frac_lo = c(0, 0.5), frac_hi = c(0.5, 1),
        density_q1 = c(0.05, 0.005), density_q2 = c(0.05, 0.005), density_q3 = c(0.05, 0.005),
        fwhm_q1 = c(5, 5), fwhm_q3 = c(10, 10)
    )
    sig <- gen_synthetic_1d_regions(n = 2000, region_profile = region_profile, seed = 1, cap_density = FALSE)
    n_dense <- sum(sig$peak_info$center <= 1000)
    n_sparse <- sum(sig$peak_info$center > 1000)
    expect_true(n_dense > n_sparse)
})

## generate_synthetic_pool_1d_regions -------------------------------------------

test_that("generate_synthetic_pool_1d_regions generates n_synthetic distinct draws", {
    region_profile <- data.frame(
        region = 1, frac_lo = 0, frac_hi = 1,
        density_q1 = 0.02, density_q2 = 0.02, density_q3 = 0.02, fwhm_q1 = 10, fwhm_q3 = 30
    )
    pooled_stats <- list(csnr = 0.03, regions = region_profile)
    pool <- generate_synthetic_pool_1d_regions(pooled_stats, list(seq_len(500)), n_synthetic = 3, peak_shape = "lorentzian")
    expect_equal(length(pool), 3)
    expect_false(identical(pool[[1]]$y, pool[[2]]$y))
})

## tune_psalsa(num_regions = ...) end-to-end ------------------------------------

test_that("tune_psalsa with num_regions set runs and returns usable parameters", {
    set.seed(3)
    x <- seq_len(400)
    baseline <- 10 + 3 * sin(x / 80)
    peaks <- 40 * exp(-((x - 60)^2) / (2 * 3^2)) + 40 * exp(-((x - 90)^2) / (2 * 3^2)) +
        15 * exp(-((x - 300)^2) / (2 * 5^2))
    y <- baseline + peaks + rnorm(length(x), 0, 0.5)

    result <- tune_psalsa(y, num_regions = 8, n_synthetic = 3, optim_maxit = 30)
    expect_true(is.finite(result$lambda) && result$lambda > 0)
    expect_true(is.finite(result$p) && result$p > 0 && result$p < 0.5)
    expect_true(is.finite(result$k) && result$k > 0)
    expect_equal(length(result$baseline), length(y))
})

test_that("tune_psalsa's default (num_regions = NULL) still uses the signal-wide pooling path", {
    set.seed(4)
    x <- seq_len(300)
    baseline <- 5 + 2 * sin(x / 60)
    peaks <- 30 * exp(-((x - 150)^2) / (2 * 4^2))
    y <- baseline + peaks + rnorm(length(x), 0, 0.3)

    result <- tune_psalsa(y, n_synthetic = 3, optim_maxit = 30)
    expect_true(is.finite(result$lambda) && result$lambda > 0)
    expect_equal(length(result$baseline), length(y))
})
