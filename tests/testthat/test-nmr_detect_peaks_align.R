make_dataset <- function(with_baseline = FALSE) {
    ds <- new_nmr_dataset_1D(
        ppm_axis = seq(0, 10, length.out = 20),
        data_1r = matrix(c(1:20, 21:40), nrow = 2, byrow = TRUE),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20")))
    )
    if (with_baseline) {
        ds$data_1r_baseline <- matrix(1, nrow = 2, ncol = 20)
    }
    ds
}

## peakList_to_dataframe / peak_data_to_peakList (round trip) -----------------

test_that("peakList_to_dataframe builds one row per peak with ppm/pos/intensity columns", {
    ds <- make_dataset()
    peakList <- list(c(3, 7), 5)

    df <- peakList_to_dataframe(ds, peakList)

    expect_equal(nrow(df), 3)
    expect_equal(df$NMRExperiment, c("10", "10", "20"))
    expect_equal(df$pos, c(3, 7, 5))
    expect_equal(df$ppm, ds$axis[c(3, 7, 5)])
    # No baseline: intensity == intensity_raw == the raw spectrum value:
    expect_equal(df$intensity_raw, c(3, 7, 25))
    expect_equal(df$intensity, df$intensity_raw)
    # peak_id is zero-padded to the total peak count's digit width:
    expect_equal(df$peak_id, c("Peak1", "Peak2", "Peak3"))
})

test_that("peakList_to_dataframe subtracts the baseline from intensity when present", {
    ds <- make_dataset(with_baseline = TRUE)
    peakList <- list(3, 5)

    df <- peakList_to_dataframe(ds, peakList)

    expect_equal(df$intensity_raw, c(3, 25))
    expect_equal(df$intensity, c(2, 24)) # baseline is a constant 1
})

test_that("peakList_to_dataframe zero-pads peak_id to the total peak count's width", {
    ds <- make_dataset()
    peakList <- list(1:12, numeric(0))

    df <- peakList_to_dataframe(ds, peakList)

    expect_equal(df$peak_id, sprintf("Peak%02d", 1:12))
})

test_that("peak_data_to_peakList inverts peakList_to_dataframe", {
    ds <- make_dataset()
    peakList <- list(c(3, 7), 5)

    df <- peakList_to_dataframe(ds, peakList)
    roundtrip <- peak_data_to_peakList(ds, df)

    expect_equal(roundtrip, peakList)
})

test_that("peak_data_to_peakList assigns an empty numeric vector for samples with no peaks", {
    ds <- make_dataset()
    peak_data <- data.frame(NMRExperiment = "10", pos = 3)

    result <- peak_data_to_peakList(ds, peak_data)

    expect_equal(result, list(3, numeric(0)))
})

## regions_from_peak_table -----------------------------------------------------

test_that("regions_from_peak_table centers a region of the given width on each peak", {
    result <- regions_from_peak_table(c(1, 2, 3), 0.1)
    expect_equal(result, list(c(0.95, 1.05), c(1.95, 2.05), c(2.95, 3.05)))
})

test_that("regions_from_peak_table accepts a per-peak width vector", {
    result <- regions_from_peak_table(c(1, 2, 3), c(0.1, 0.2, 0.3))
    expect_equal(result, list(c(0.95, 1.05), c(1.9, 2.1), c(2.85, 3.15)))
})

## nmr_ppm_resolution / ppm_resolution ----------------------------------------

test_that("nmr_ppm_resolution.nmr_dataset_1D returns the median spacing of the ppm axis", {
    ds <- new_nmr_dataset_1D(
        ppm_axis = seq(0, 10, length.out = 21), # spacing = 0.5
        data_1r = matrix(1:42, nrow = 2, byrow = TRUE),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20")))
    )

    expect_equal(nmr_ppm_resolution(ds), 0.5)
})

test_that("ppm_resolution unlists nmr_ppm_resolution() on the first sample", {
    ds <- new_nmr_dataset_1D(
        ppm_axis = seq(0, 10, length.out = 21),
        data_1r = matrix(1:42, nrow = 2, byrow = TRUE),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20")))
    )

    expect_equal(ppm_resolution(ds), 0.5)
    expect_type(ppm_resolution(ds), "double")
})

## signif_transformer (glue transformer factory) -------------------------------

test_that("signif_transformer rounds numeric glue substitutions to the given number of significant digits", {
    tr <- signif_transformer(2)
    val <- 3.14159

    result <- glue::glue("{val}", .transformer = tr)

    expect_equal(as.character(result), "3.1")
})

test_that("signif_transformer leaves non-numeric glue substitutions unchanged", {
    tr <- signif_transformer(2)
    name <- "sample"

    result <- glue::glue("{name}", .transformer = tr)

    expect_equal(as.character(result), "sample")
})

## nmr_detect_peaks_plot_overview -----------------------------------------------

test_that("nmr_detect_peaks_plot_overview filters to accepted peaks by default", {
    skip_if_not_installed("ggplot2")
    peak_data <- data.frame(
        NMRExperiment = c("10", "10", "20"),
        ppm = c(1.52, 3.18, 1.63),
        accepted = c(TRUE, TRUE, FALSE)
    )

    gplt <- nmr_detect_peaks_plot_overview(peak_data)

    expect_s3_class(gplt, "ggplot")
    expect_equal(sum(gplt$data$num_peaks), 2)
})

test_that("nmr_detect_peaks_plot_overview can include non-accepted peaks", {
    skip_if_not_installed("ggplot2")
    peak_data <- data.frame(
        NMRExperiment = c("10", "10", "20"),
        ppm = c(1.52, 3.18, 1.63),
        accepted = c(TRUE, TRUE, FALSE)
    )

    gplt <- nmr_detect_peaks_plot_overview(peak_data, accepted_only = FALSE)

    expect_equal(sum(gplt$data$num_peaks), 3)
})
