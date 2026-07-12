make_dataset <- function(n = 2, seed = 1, with_baseline = FALSE) {
    set.seed(seed)
    ppm_axis <- seq(from = 0, to = 10, length.out = 1000)
    data_1r <- matrix(rnorm(n * 1000, mean = 100, sd = 2), nrow = n)
    dataset_1D <- new_nmr_dataset_1D(
        ppm_axis = ppm_axis,
        data_1r = data_1r,
        metadata = list(external = data.frame(NMRExperiment = as.character(seq_len(n) * 10)))
    )
    if (with_baseline) {
        dataset_1D$data_1r_baseline <- matrix(5, nrow = n, ncol = 1000)
    }
    dataset_1D
}

## nmr_baseline_threshold ------------------------------------------------------

test_that("nmr_baseline_threshold's mean3sd method matches a manual mean+3sd calculation (multiple samples)", {
    dataset_1D <- make_dataset(n = 2)

    result <- nmr_baseline_threshold(dataset_1D, range_without_peaks = c(9.5, 10), method = "mean3sd")

    threshold_ind <- dataset_1D$axis >= 9.5 & dataset_1D$axis < 10
    region <- dataset_1D$data_1r[, threshold_ind, drop = FALSE]
    expected <- mean(apply(region, 2, mean)) + 3 * mean(apply(region, 2, stats::sd))
    expect_equal(result, expected)
    # mean3sd returns a single, unnamed scalar for the whole dataset:
    expect_length(result, 1)
})

test_that("nmr_baseline_threshold's mean3sd method matches a manual mean+3sd calculation (single sample)", {
    dataset_1D <- make_dataset(n = 1)

    result <- nmr_baseline_threshold(dataset_1D, range_without_peaks = c(9.5, 10), method = "mean3sd")

    threshold_ind <- dataset_1D$axis >= 9.5 & dataset_1D$axis < 10
    region <- as.numeric(dataset_1D$data_1r[, threshold_ind])
    expected <- mean(region) + 3 * stats::sd(region)
    expect_equal(result, expected)
})

test_that("nmr_baseline_threshold's median3mad method returns one named value per sample", {
    dataset_1D <- make_dataset(n = 2)

    result <- nmr_baseline_threshold(dataset_1D, range_without_peaks = c(9.5, 10), method = "median3mad")

    threshold_ind <- dataset_1D$axis >= 9.5 & dataset_1D$axis < 10
    expected <- vapply(seq_len(2), function(i) {
        region <- dataset_1D$data_1r[i, threshold_ind, drop = FALSE]
        stats::median(region) + 3 * stats::mad(region)
    }, numeric(1))
    expect_equal(unname(result), expected)
    expect_equal(names(result), names(dataset_1D))
})

test_that("nmr_baseline_threshold's median3mad method subtracts data_1r_baseline when present", {
    dataset_1D <- make_dataset(n = 2, with_baseline = TRUE)

    result <- nmr_baseline_threshold(dataset_1D, range_without_peaks = c(9.5, 10), method = "median3mad")

    threshold_ind <- dataset_1D$axis >= 9.5 & dataset_1D$axis < 10
    expected <- vapply(seq_len(2), function(i) {
        region <- dataset_1D$data_1r[i, threshold_ind, drop = FALSE] - dataset_1D$data_1r_baseline[i, threshold_ind]
        stats::median(region) + 3 * stats::mad(region)
    }, numeric(1))
    expect_equal(unname(result), expected)

    # The baseline-subtracted threshold must differ from the raw one (the
    # synthetic baseline here is a nonzero constant offset):
    result_no_baseline <- nmr_baseline_threshold(
        make_dataset(n = 2, with_baseline = FALSE),
        range_without_peaks = c(9.5, 10),
        method = "median3mad"
    )
    expect_false(isTRUE(all.equal(unname(result), unname(result_no_baseline))))
})

test_that("nmr_baseline_threshold requires range_without_peaks to be given", {
    dataset_1D <- make_dataset(n = 1)
    expect_error(
        nmr_baseline_threshold(dataset_1D),
        "range_without_peaks must be given"
    )
})

test_that("nmr_baseline_threshold requires range_without_peaks to have length 2", {
    dataset_1D <- make_dataset(n = 1)
    expect_error(
        nmr_baseline_threshold(dataset_1D, range_without_peaks = c(1, 2, 3)),
        "length 2"
    )
})

test_that("nmr_baseline_threshold refuses to estimate a threshold from too few points", {
    dataset_1D <- make_dataset(n = 1)
    expect_error(
        nmr_baseline_threshold(dataset_1D, range_without_peaks = c(9.999, 10)),
        "Can't estimate a baseline threshold reliably"
    )
})

## tidy_spectra_baseline_and_threshold (internal helper) ----------------------

test_that("tidy_spectra_baseline_and_threshold returns NULL baselines/thresholds when neither is available", {
    dataset_1D <- make_dataset(n = 2)

    result <- tidy_spectra_baseline_and_threshold(
        dataset_1D,
        thresholds = NULL,
        chemshift_range = c(9.5, 10),
        NMRExperiment = c("10", "20")
    )

    expect_named(result, c("spectra", "baselines", "thresholds"))
    expect_null(result$baselines)
    expect_null(result$thresholds)
    expect_true(all(c("NMRExperiment", "chemshift", "intensity") %in% colnames(result$spectra)))
})

test_that("tidy_spectra_baseline_and_threshold reports the raw threshold value when there is no baseline", {
    dataset_1D <- make_dataset(n = 2)
    th <- nmr_baseline_threshold(dataset_1D, range_without_peaks = c(9.5, 10), method = "median3mad")

    result <- tidy_spectra_baseline_and_threshold(
        dataset_1D,
        thresholds = th,
        chemshift_range = c(9.5, 10),
        NMRExperiment = c("10", "20")
    )

    expect_null(result$baselines)
    expect_true(all(result$thresholds$intensity == th[result$thresholds$NMRExperiment]))
})

test_that("tidy_spectra_baseline_and_threshold offsets the threshold by the baseline when one is present", {
    dataset_1D <- make_dataset(n = 2, with_baseline = TRUE)
    th <- nmr_baseline_threshold(dataset_1D, range_without_peaks = c(9.5, 10), method = "median3mad")

    result <- tidy_spectra_baseline_and_threshold(
        dataset_1D,
        thresholds = th,
        chemshift_range = c(9.5, 10),
        NMRExperiment = c("10", "20")
    )

    expect_equal(unique(result$baselines$intensity), 5)
    # threshold trace = baseline (5) + the per-sample threshold value:
    expect_true(all(result$thresholds$intensity == 5 + th[result$thresholds$NMRExperiment]))
})

## nmr_baseline_threshold_plot --------------------------------------------------

test_that("nmr_baseline_threshold_plot requires chemshift_range to be given", {
    dataset_1D <- make_dataset(n = 1)
    th <- nmr_baseline_threshold(dataset_1D, range_without_peaks = c(9.5, 10), method = "median3mad")
    expect_error(
        nmr_baseline_threshold_plot(dataset_1D, th),
        "chemshift_range must be given"
    )
})

test_that("nmr_baseline_threshold_plot returns a ggplot for the default NMRExperiment='all'", {
    skip_if_not_installed("ggplot2")
    dataset_1D <- make_dataset(n = 2)
    th <- nmr_baseline_threshold(dataset_1D, range_without_peaks = c(9.5, 10), method = "median3mad")

    gplt <- nmr_baseline_threshold_plot(dataset_1D, th, NMRExperiment = "all", chemshift_range = c(9.5, 10))

    expect_s3_class(gplt, "ggplot")
})

test_that("nmr_baseline_threshold_plot recycles a single threshold value across all requested experiments", {
    skip_if_not_installed("ggplot2")
    dataset_1D <- make_dataset(n = 2)

    gplt <- nmr_baseline_threshold_plot(dataset_1D, thresholds = 105, NMRExperiment = "all", chemshift_range = c(9.5, 10))

    expect_s3_class(gplt, "ggplot")
})

test_that("nmr_baseline_threshold_plot subsets thresholds to the requested NMRExperiment", {
    skip_if_not_installed("ggplot2")
    dataset_1D <- make_dataset(n = 2)
    th <- nmr_baseline_threshold(dataset_1D, range_without_peaks = c(9.5, 10), method = "median3mad")

    gplt <- nmr_baseline_threshold_plot(dataset_1D, th, NMRExperiment = "10", chemshift_range = c(9.5, 10))

    expect_s3_class(gplt, "ggplot")
    expect_equal(unique(gplt$layers[[1]]$data$NMRExperiment), "10")
})

test_that("nmr_baseline_threshold_plot warns once when passed deprecated aes_string arguments", {
    skip_if_not_installed("ggplot2")
    dataset_1D <- make_dataset(n = 2)
    th <- nmr_baseline_threshold(dataset_1D, range_without_peaks = c(9.5, 10), method = "median3mad")

    # The deprecation warning is throttled to once per cli .frequency_id per
    # session (see nmr_baseline_threshold_plot()'s use of .frequency =
    # "regularly"), so only exercise this path in a single test.
    expect_warning(
        gplt <- nmr_baseline_threshold_plot(dataset_1D, th, NMRExperiment = "all", chemshift_range = c(9.5, 10), colour = "NMRExperiment"),
        "deprecated"
    )
    expect_s3_class(gplt, "ggplot")
})
