test_that("nmr_normalize & nmr_normalize_extra_info work", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:10,
        data_1r = matrix(sample(0:99, replace = TRUE), nrow = 10),
        metadata = list(external = data.frame(
            NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100")
        ))
    )
    dataset <- nmr_normalize(dataset, method = "pqn")
    diagnostic <- nmr_normalize_extra_info(dataset)
    expect_true(is.matrix(dataset[["data_1r"]]))
    expect_true(is.numeric(dataset[["data_1r"]][[1]]))
    expect_true(is.list(diagnostic))
    expect_true(is.data.frame(diagnostic[["norm_factor"]]))
    expect_true(is.character(diagnostic[["norm_factor"]][[1, 1]]))
})

test_that("nmr_normalize with negative values and below 10 samples work", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:10,
        data_1r = matrix(sample(-1:88, replace = FALSE), nrow = 9),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90")))
    )
    expect_warning(nmr_normalize(dataset, method = "pqn"))
})

test_that("nmr_normalize works with unknown method", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:10,
        data_1r = matrix(sample(-1:88, replace = FALSE), nrow = 9),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90")))
    )
    dataset_none <- nmr_normalize(dataset, method = "none")
    dataset_area <- nmr_normalize(dataset, method = "area")
    dataset_max <- nmr_normalize(dataset, method = "max")
    dataset_ <- nmr_normalize(dataset)

    expect_true(is.matrix(dataset_none[["data_1r"]]))
    expect_true(is.matrix(dataset_area[["data_1r"]]))
    expect_true(is.matrix(dataset_max[["data_1r"]]))
    expect_true(is.matrix(dataset_[["data_1r"]]))
})

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
    # isolates the median == 0 check.
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

test_that("nmr_normalize warns when a normalization factor is non-positive", {
    # apply(., 1, max): 10, 20, -1 -> one non-positive factor, but
    # median(c(10, 20, -1)) == 10, so the median == 0 check does NOT fire.
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
