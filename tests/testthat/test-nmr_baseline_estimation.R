## Small, fast synthetic dataset shared across these tests: 3 samples, each a
## smooth baseline plus a couple of peaks plus a little noise.
make_test_dataset <- function(seed = 1, n = 300, nsamples = 3) {
    set.seed(seed)
    x <- seq_len(n)
    ppm_axis <- seq(0.5, 9.5, length.out = n)
    spectra <- t(vapply(seq_len(nsamples), function(i) {
        baseline <- 5 + 2 * sin(x / 80)
        peaks <- 20 * exp(-((x - 100)^2) / (2 * 4^2)) + 15 * exp(-((x - 220)^2) / (2 * 5^2))
        baseline + peaks + stats::rnorm(n, 0, 0.3)
    }, numeric(n)))
    new_nmr_dataset_1D(
        ppm_axis = ppm_axis,
        data_1r = spectra,
        metadata = list(external = data.frame(NMRExperiment = paste0("Sample", seq_len(nsamples))))
    )
}

test_that("nmr_baseline_estimation defaults to scalar tune_psalsa (num_regions = NULL)", {
    dataset <- make_test_dataset()
    result <- nmr_baseline_estimation(dataset)
    expect_true("data_1r_baseline" %in% names(unclass(result)))
    expect_equal(dim(result$data_1r_baseline), dim(dataset$data_1r))

    params <- attr(result$data_1r_baseline, "psalsa_params")
    expect_equal(length(params$lambda), 1)
    expect_equal(length(params$p), 1)
    expect_equal(length(params$k), 1)
    expect_null(params$num_regions)
})

## A denser fixture, needed specifically for spatial tuning to actually
## produce a multi-region profile instead of falling back to a scalar tune:
## tune_psalsa_spatial()'s merging (see merge_sparse_regions_1d()) collapses
## regions with too few real peaks into one another, and with the default
## min_peaks = 15, a handful of peaks in a 300-point toy signal (as in
## make_test_dataset()) isn't enough to keep more than one group. This one
## packs >= 15 detectable peaks into each of 2 of its 4 quarters.
make_dense_test_dataset <- function(seed = 1, n = 800, nsamples = 3) {
    set.seed(seed)
    x <- seq_len(n)
    ppm_axis <- seq(0.5, 9.5, length.out = n)
    peak_centers <- c(seq(10, 190, by = 10), seq(210, 390, by = 10))
    spectra <- t(vapply(seq_len(nsamples), function(i) {
        baseline <- 5 + 2 * sin(x / 200)
        peaks <- Reduce(`+`, lapply(peak_centers, function(ctr) {
            15 * exp(-((x - ctr)^2) / (2 * 1.5^2))
        }))
        baseline + peaks + stats::rnorm(n, 0, 0.2)
    }, numeric(n)))
    new_nmr_dataset_1D(
        ppm_axis = ppm_axis,
        data_1r = spectra,
        metadata = list(external = data.frame(NMRExperiment = paste0("Sample", seq_len(nsamples))))
    )
}

test_that("nmr_baseline_estimation with num_regions set tunes a position-varying profile", {
    dataset <- make_dense_test_dataset()
    result <- nmr_baseline_estimation(dataset, num_regions = 4)
    expect_true("data_1r_baseline" %in% names(unclass(result)))
    expect_equal(dim(result$data_1r_baseline), dim(dataset$data_1r))

    params <- attr(result$data_1r_baseline, "psalsa_params")
    n <- ncol(dataset$data_1r)
    expect_equal(length(params$lambda), n)
    expect_equal(length(params$p), n)
    expect_equal(length(params$k), n)
    expect_equal(params$num_regions, 4)
    expect_true(all(is.finite(result$data_1r_baseline)))
})

test_that("explicit numeric lambda/p/k bypass tuning even when num_regions is set", {
    dataset <- make_test_dataset()
    result <- nmr_baseline_estimation(dataset, lambda = 1e6, p = 0.01, k = 5, num_regions = 10)
    params <- attr(result$data_1r_baseline, "psalsa_params")
    expect_equal(params$lambda, 1e6)
    expect_equal(params$p, 0.01)
    expect_equal(params$k, 5)
    # num_regions is still recorded as given, even though it had no effect here:
    expect_equal(params$num_regions, 10)
})

test_that("num_regions is ignored when lambda/p/k are all given explicitly", {
    dataset <- make_test_dataset()
    # num_regions = 4 would normally trigger tune_psalsa_spatial(), but there's
    # nothing left to tune -- this must not error or attempt any tuning.
    result <- nmr_baseline_estimation(dataset, lambda = 1e6, p = 0.01, k = 5, num_regions = 4)
    expect_true(all(is.finite(result$data_1r_baseline)))
})

test_that("psalsa_params from a num_regions tune can be reused directly on a same-sized dataset", {
    dataset1 <- make_dense_test_dataset(seed = 10)
    tuned <- nmr_baseline_estimation(dataset1, num_regions = 4)
    params <- attr(tuned$data_1r_baseline, "psalsa_params")

    dataset2 <- make_dense_test_dataset(seed = 11)
    reused <- nmr_baseline_estimation(dataset2, lambda = params$lambda, p = params$p, k = params$k)
    expect_equal(dim(reused$data_1r_baseline), dim(dataset2$data_1r))
    expect_true(all(is.finite(reused$data_1r_baseline)))
    # No tuning happened on reuse -- the params attribute is exactly what was passed in:
    expect_equal(attr(reused$data_1r_baseline, "psalsa_params")$lambda, params$lambda)
})
