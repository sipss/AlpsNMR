## peak_is_local_max_1 -----------------------------------------------------

test_that("peak_is_local_max_1 accepts a genuine local maximum", {
    # values: 1 2 3 4 5 4 3 2 1, apex at index 5 (value 5), window [3,7]
    intensity <- c(1, 2, 3, 4, 5, 4, 3, 2, 1)
    expect_true(peak_is_local_max_1(intensity, pos = 5, idx_lo = 3, idx_hi = 7, tol_frac = 0.01))
})

test_that("peak_is_local_max_1 rejects a point on a monotonic ramp (shoulder)", {
    # a plain monotonic ramp, as if riding the flank of a much wider,
    # neighbouring peak whose true apex sits outside this window
    intensity <- 1:11
    expect_false(peak_is_local_max_1(intensity, pos = 6, idx_lo = 4, idx_hi = 8, tol_frac = 0.01))
})

test_that("peak_is_local_max_1 accepts a near-tie within tol_frac (flat/rounded apex)", {
    # true max (30) is one point to the right of pos, but well within a 1% tolerance
    intensity <- c(10, 20, 29.9, 30, 29.8, 20, 10)
    expect_true(peak_is_local_max_1(intensity, pos = 3, idx_lo = 1, idx_hi = 7, tol_frac = 0.01))
})

test_that("peak_is_local_max_1 rejects a gap larger than tol_frac", {
    # pos is 10% below the window's true max -- well past a 1% tolerance
    intensity <- c(10, 20, 27, 30, 27, 20, 10)
    expect_false(peak_is_local_max_1(intensity, pos = 3, idx_lo = 1, idx_hi = 7, tol_frac = 0.01))
})

test_that("peak_is_local_max_1 handles an entirely negative window (regression test)", {
    # A window whose intensity dips below zero (e.g. baseline noise). A
    # tolerance computed as tol_frac * wmax without abs() would go negative
    # here and invert every comparison -- even an exact tie would then be
    # rejected. wmax = -8 at index 4 (pos), an exact match.
    intensity <- c(-10, -9, -8.1, -8, -8.05, -8.3, -8.9, -9.5, -10.2)
    expect_true(peak_is_local_max_1(intensity, pos = 4, idx_lo = 1, idx_hi = 9, tol_frac = 0.01))
})

test_that("peak_is_local_max_1 returns NA for a degenerate window", {
    intensity <- 1:11
    # window too small (fewer than 3 points)
    expect_true(is.na(peak_is_local_max_1(intensity, pos = 5, idx_lo = 5, idx_hi = 6, tol_frac = 0.01)))
    # pos outside the window
    expect_true(is.na(peak_is_local_max_1(intensity, pos = 2, idx_lo = 4, idx_hi = 8, tol_frac = 0.01)))
})

## is_peak_local_max ---------------------------------------------------------

test_that("is_peak_local_max classifies peaks per-sample using each sample's own raw intensity", {
    nmr_dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:11,
        data_1r = rbind(
            c(1, 2, 3, 4, 5, 4, 3, 2, 1, 1, 1), # sample "10": genuine peak at ppm 5
            c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11) # sample "20": monotonic ramp
        ),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20")))
    )
    peak_data <- data.frame(
        peak_id = c("PeakA", "PeakB"),
        NMRExperiment = c("10", "20"),
        pos = c(5, 6),
        ppm_infl_min = c(3, 4),
        ppm_infl_max = c(7, 8)
    )
    result <- is_peak_local_max(peak_data, nmr_dataset)
    expect_equal(result, c(TRUE, FALSE))
})

## peaklist_accept_peaks: accept_inflections ----------------------------------

test_that("peaklist_accept_peaks's accept_inflections defaults to TRUE (previous behaviour)", {
    nmr_dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:11,
        data_1r = matrix(1:11, nrow = 1),
        metadata = list(external = data.frame(NMRExperiment = "10"))
    )
    peak_data <- data.frame(
        peak_id = "PeakA",
        NMRExperiment = "10",
        ppm = 6,
        pos = 6,
        intensity = 6,
        ppm_infl_min = 4,
        ppm_infl_max = 8,
        area = 100,
        norm_rmse = 0.01
    )
    # A shoulder on a monotonic ramp, but accept_inflections isn't set, so
    # it's accepted like any other peak that passes the other criteria.
    result <- peaklist_accept_peaks(peak_data, nmr_dataset, area_min = 10)
    expect_true(result$accepted)
})

test_that("peaklist_accept_peaks with accept_inflections = FALSE rejects shoulders but keeps real peaks", {
    nmr_dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:11,
        data_1r = rbind(
            c(1, 2, 3, 4, 5, 4, 3, 2, 1, 1, 1),
            c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
        ),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20")))
    )
    peak_data <- data.frame(
        peak_id = c("PeakA", "PeakB"),
        NMRExperiment = c("10", "20"),
        ppm = c(5, 6),
        pos = c(5, 6),
        intensity = c(5, 6),
        ppm_infl_min = c(3, 4),
        ppm_infl_max = c(7, 8),
        area = c(100, 100),
        norm_rmse = c(0.01, 0.01)
    )
    result <- peaklist_accept_peaks(peak_data, nmr_dataset, area_min = 10, accept_inflections = FALSE)
    expect_equal(result$accepted, c(TRUE, FALSE))
})

test_that("peaklist_accept_peaks with accept_inflections = FALSE and keep_rejected = FALSE drops shoulders", {
    nmr_dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:11,
        data_1r = rbind(
            c(1, 2, 3, 4, 5, 4, 3, 2, 1, 1, 1),
            c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
        ),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20")))
    )
    peak_data <- data.frame(
        peak_id = c("PeakA", "PeakB"),
        NMRExperiment = c("10", "20"),
        ppm = c(5, 6),
        pos = c(5, 6),
        intensity = c(5, 6),
        ppm_infl_min = c(3, 4),
        ppm_infl_max = c(7, 8),
        area = c(100, 100),
        norm_rmse = c(0.01, 0.01)
    )
    result <- peaklist_accept_peaks(
        peak_data, nmr_dataset,
        area_min = 10, accept_inflections = FALSE, keep_rejected = FALSE
    )
    expect_equal(result$peak_id, "PeakA")
    expect_false("accepted" %in% colnames(result))
})

## peaklist_accept_peaks: intensity_min / intensity_max -----------------------

test_that("peaklist_accept_peaks's intensity_min/intensity_max reject on height, not area", {
    nmr_dataset <- new_nmr_dataset_1D(
        1:10,
        matrix(c(1:5, 4:2, 3, 0), nrow = 1),
        list(external = data.frame(NMRExperiment = "10"))
    )
    peak_data <- data.frame(
        peak_id = c("Tall", "Short"),
        NMRExperiment = c("10", "10"),
        ppm = c(5, 9),
        pos = c(5, 9),
        intensity = c(100, 3),
        ppm_infl_min = c(3, 8),
        ppm_infl_max = c(7, 10),
        gamma_ppb = c(1, 1),
        # Areas deliberately don't track intensity, to confirm the new
        # criterion filters on intensity and not area:
        area = c(3, 100),
        norm_rmse = c(0.01, 0.01)
    )

    by_intensity_min <- peaklist_accept_peaks(peak_data, nmr_dataset, intensity_min = 10)
    expect_equal(by_intensity_min$accepted, c(TRUE, FALSE))

    by_intensity_max <- peaklist_accept_peaks(peak_data, nmr_dataset, intensity_max = 10)
    expect_equal(by_intensity_max$accepted, c(FALSE, TRUE))
})

test_that("peaklist_accept_peaks defaults intensity_min/intensity_max to a no-op", {
    nmr_dataset <- new_nmr_dataset_1D(
        1:10,
        matrix(c(1:5, 4:2, 3, 0), nrow = 1),
        list(external = data.frame(NMRExperiment = "10"))
    )
    peak_data <- data.frame(
        peak_id = c("Peak1", "Peak2"),
        NMRExperiment = c("10", "10"),
        ppm = c(5, 9),
        pos = c(5, 9),
        intensity = c(1e6, 1e-6),
        ppm_infl_min = c(3, 8),
        ppm_infl_max = c(7, 10),
        gamma_ppb = c(1, 1),
        area = c(25, 25),
        norm_rmse = c(0.01, 0.01)
    )
    result <- peaklist_accept_peaks(peak_data, nmr_dataset)
    expect_equal(result$accepted, c(TRUE, TRUE))
})
