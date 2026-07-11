# Regression tests for a set of bug fixes identified in a critical code
# review of R/nmr_autophase.R, R/nmr_normalize.R and R/nmr_interpolate.R.
#
# Each test_that() block below is annotated with the specific fix it
# exercises so that a regression in any of these fixed behaviors is caught.

# ---------------------------------------------------------------------------
# nmr_autophase(): absorptionOnly passthrough fix
#
# Previously the call `NMRphasing::NMRphasing(to_phase, absorptionOnly = TRUE, ...)`
# hardcoded `TRUE`, discarding the imaginary/phase information even when it
# was available. It is now `absorptionOnly = absorptionOnly`, the
# locally-computed per-sample variable (FALSE when imaginary data is present,
# TRUE otherwise).
#
# We verify this end-to-end by mocking NMRphasing::NMRphasing() and recording
# the actual value of `absorptionOnly` it was called with, for a dataset
# where one sample has imaginary data and another does not.
# ---------------------------------------------------------------------------
test_that("nmr_autophase forwards the locally-computed absorptionOnly to NMRphasing::NMRphasing", {
    skip_if_not_installed("NMRphasing")

    # Force serial execution so the mocked binding (only visible in this
    # process) is actually used, and so results are recorded in a
    # deterministic, per-sample order.
    old_bpparam <- BiocParallel::bpparam()
    BiocParallel::register(old_bpparam)
    BiocParallel::register(BiocParallel::SerialParam())
    on.exit(BiocParallel::register(old_bpparam), add = TRUE)

    log_file <- tempfile()
    on.exit(unlink(log_file), add = TRUE)
    testthat::local_mocked_bindings(
        NMRphasing = function(X, absorptionOnly, ...) {
            # Record what nmr_autophase() actually asked for.
            cat(absorptionOnly, "\n", file = log_file, append = TRUE)
            X
        },
        .package = "NMRphasing"
    )

    x <- seq(from = 1, to = 2, length.out = 20)
    y_real <- x
    y_imag <- x * 2
    dataset <- new_nmr_dataset(
        metadata = list(external = data.frame(NMRExperiment = c("10", "20"))),
        data_fields = list(
            # Sample "10" has imaginary data, sample "20" does not.
            data_1r = list(y_real, y_real),
            data_1i = list(y_imag, NULL)
        ),
        axis = list(list(x), list(x))
    )

    expect_warning(
        nmr_autophase(dataset, method = "NLS"),
        "missing imaginary spectrum"
    )

    logged <- scan(log_file, what = character(), quiet = TRUE)
    expect_length(logged, 2)
    # Sample with imaginary data => absorptionOnly must be FALSE (full
    # complex phasing). Sample without it => absorptionOnly must be TRUE.
    # Before the fix, both calls would have received TRUE.
    expect_equal(logged, c("FALSE", "TRUE"))
})

# ---------------------------------------------------------------------------
# nmr_autophase(): misplaced-parenthesis fix in the "missing imaginary
# component" message.
#
# Previously `if (length(any_imag_missing < 7))` compared a logical vector
# elementwise to 7 (always a positive-length, always-truthy result), which
# made the "truncate to first 5 + N more" branch dead code: every run with
# missing imaginary data listed every sample name, however many there were.
# It is now `if (length(any_imag_missing) < 7)`, so >= 7 missing samples
# trigger the "first 5 ... and N more" truncated form.
#
# NMRphasing::NMRphasing() is mocked here (as an identity function) purely
# to keep the test fast and independent of the phasing algorithm's numerical
# behavior on degenerate synthetic data -- we only care about the message
# that is built *before* NMRphasing is invoked.
# ---------------------------------------------------------------------------
test_that("nmr_autophase summarizes many missing-imaginary samples instead of listing them all", {
    skip_if_not_installed("NMRphasing")

    testthat::local_mocked_bindings(
        NMRphasing = function(X, absorptionOnly, ...) X,
        .package = "NMRphasing"
    )

    old_width <- options(width = 1000)
    on.exit(options(old_width), add = TRUE)

    n <- 10
    x <- seq(from = 1, to = 2, length.out = 20)
    y <- x
    sample_names <- paste0("Sample", seq_len(n))
    dataset <- new_nmr_dataset(
        metadata = list(external = data.frame(NMRExperiment = sample_names)),
        data_fields = list(
            data_1r = stats::setNames(rep(list(y), n), sample_names),
            # None of the 10 samples has imaginary data.
            data_1i = stats::setNames(vector("list", n), sample_names)
        ),
        axis = rep(list(list(x)), n)
    )

    warnings_raised <- testthat::capture_warnings(
        nmr_autophase(dataset, method = "NLS")
    )
    combined <- paste(warnings_raised, collapse = " | ")

    # The truncated form must be used ...
    expect_match(combined, "and 5 more", fixed = TRUE)
    # ... and not the full dump of all 10 sample names (in particular the
    # 6th-to-10th names must not appear).
    for (missing_name in sample_names[6:10]) {
        expect_false(grepl(missing_name, combined, fixed = TRUE))
    }
})

# ---------------------------------------------------------------------------
# nmr_normalize() / norm_pqn(): divide-by-zero guard.
#
# Previously `norm_factor_to_apply <- norm_factor / stats::median(norm_factor)`
# had no guard against `median(norm_factor) == 0`, silently producing
# Inf/NaN spectra. There is now an explicit check that calls stop() with a
# clear message.
# ---------------------------------------------------------------------------
test_that("nmr_normalize errors clearly when the median normalization factor is zero (method = 'area')", {
    # Row sums (the "area" norm_factor): -25, -5, 5, 25 -> median is exactly 0.
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:5,
        data_1r = rbind(
            c(-5, -5, -5, -5, -5),
            c(-1, -1, -1, -1, -1),
            c(1, 1, 1, 1, 1),
            c(5, 5, 5, 5, 5)
        ),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30", "40")))
    )

    # Some norm_factor values are also <= 0 here, which independently raises
    # the non-positive-factor warning tested below; suppress it so this test
    # isolates the median == 0 guard.
    expect_error(
        suppressWarnings(nmr_normalize(dataset, method = "area")),
        regexp = "median of the normalization factors is zero"
    )
})

test_that("norm_pqn() errors clearly when the median normalization area is zero", {
    # rowSums: -25, -5, 5, 25 -> median area is exactly 0.
    spectra <- rbind(
        c(-5, -5, -5, -5, -5),
        c(-1, -1, -1, -1, -1),
        c(1, 1, 1, 1, 1),
        c(5, 5, 5, 5, 5)
    )
    expect_error(
        suppressWarnings(norm_pqn(spectra)),
        regexp = "median of the normalization areas is zero"
    )
})

# ---------------------------------------------------------------------------
# nmr_normalize(): non-positive-factor warning.
#
# After computing norm_factor, there is now a check
# `if (any(norm_factor <= 0)) rlang::warn(...)` that warns when a
# normalization factor is non-positive (which can flip a sample's spectrum
# sign). This must fire without necessarily hitting the median == 0 stop().
# ---------------------------------------------------------------------------
test_that("nmr_normalize warns when a normalization factor is non-positive", {
    # apply(., 1, max): 10, 20, -1 -> one non-positive factor, but
    # median(c(10, 20, -1)) == 10, so the median == 0 guard does NOT fire.
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:5,
        data_1r = rbind(
            c(1, 2, 3, 4, 10),
            c(2, 4, 6, 8, 20),
            c(-5, -4, -3, -2, -1)
        ),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30")))
    )

    expect_warning(
        nmr_normalize(dataset, method = "max"),
        regexp = "non-positive"
    )
    # And it should complete normally (no error) once the warning is heeded.
    result <- suppressWarnings(nmr_normalize(dataset, method = "max"))
    expect_true(is.matrix(result[["data_1r"]]))
})

# ---------------------------------------------------------------------------
# nmr_normalize(): method = "value" validation.
#
# Previously `dots[["values"]]` was used directly with no validation. Now
# there are explicit checks for a missing `values` argument and for a
# length mismatch against the number of samples, each with a clear stop().
# ---------------------------------------------------------------------------
test_that("nmr_normalize errors clearly when method = 'value' is used without a 'values' argument", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:5,
        data_1r = rbind(c(1, 2, 3, 4, 5), c(2, 3, 4, 5, 6)),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20")))
    )
    expect_error(
        nmr_normalize(dataset, method = "value"),
        regexp = "'values' argument is required"
    )
})

test_that("nmr_normalize errors clearly when method = 'value' values length mismatches the sample count", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:5,
        data_1r = rbind(c(1, 2, 3, 4, 5), c(2, 3, 4, 5, 6), c(3, 4, 5, 6, 7)),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30")))
    )
    # 3 samples, but only 2 values provided.
    expect_error(
        nmr_normalize(dataset, method = "value", values = c(1, 2)),
        regexp = "must have length equal to the number of samples"
    )
})

test_that("nmr_normalize works correctly when method = 'value' has properly sized values", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:5,
        data_1r = rbind(c(1, 2, 3, 4, 5), c(2, 3, 4, 5, 6), c(3, 4, 5, 6, 7)),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30")))
    )
    result <- nmr_normalize(dataset, method = "value", values = c(1, 2, 3))
    expect_true(is.matrix(result[["data_1r"]]))
})

# ---------------------------------------------------------------------------
# nmr_interpolate_1D(): out-of-range warning.
#
# Previously, if the requested interpolation axis extended beyond a
# sample's native ppm range, signal::interp1(method = "spline") silently
# extrapolated with no warning. There is now an explicit check comparing
# range(list_of_ppms[[i]]) vs range(ppm_axis) that calls rlang::warn() when
# the requested range exceeds a sample's native range.
# ---------------------------------------------------------------------------
test_that("nmr_interpolate_1D warns when the requested axis exceeds a sample's native ppm range", {
    x1 <- seq(from = 1, to = 5, length.out = 50) # wide native range
    x2 <- seq(from = 2, to = 4, length.out = 50) # narrow native range
    y1 <- sin(x1)
    y2 <- sin(x2)
    dataset <- new_nmr_dataset(
        metadata = list(external = data.frame(NMRExperiment = c("10", "20"))),
        data_fields = list(data_1r = list(y1, y2)),
        axis = list(list(x1), list(x2))
    )

    expect_warning(
        nmr_interpolate_1D(dataset, axis = c(min = 1, max = 5, by = 0.1)),
        regexp = "exceeds the native ppm range"
    )
})

test_that("nmr_interpolate_1D does not warn when the requested axis is within every sample's native range", {
    x1 <- seq(from = 1, to = 5, length.out = 50)
    x2 <- seq(from = 2, to = 4, length.out = 50)
    y1 <- sin(x1)
    y2 <- sin(x2)
    dataset <- new_nmr_dataset(
        metadata = list(external = data.frame(NMRExperiment = c("10", "20"))),
        data_fields = list(data_1r = list(y1, y2)),
        axis = list(list(x1), list(x2))
    )

    # Requested range [2.5, 3.5] is inside both samples' native ranges.
    expect_no_warning(
        nmr_interpolate_1D(dataset, axis = c(min = 2.5, max = 3.5, by = 0.1))
    )
})
