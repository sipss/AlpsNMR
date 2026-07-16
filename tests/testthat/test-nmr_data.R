test_that("new_nmr_dataset_1D", {
    ds <- new_nmr_dataset_1D(
        ppm_axis = c(1, 2, 3),
        data_1r = matrix(c(1, 2, 1), nrow = 1),
        metadata = list(external = data.frame(NMRExperiment = "10"))
    )
    dataset <- nmr_data(ds)
    expect_true(is.matrix(dataset))
})

test_that("nmr_data.nmr_dataset_1D names the rows/columns from the sample names and ppm axis", {
    ds <- new_nmr_dataset_1D(
        ppm_axis = c(1, 2, 3),
        data_1r = matrix(c(1, 2, 1, 4, 5, 6), nrow = 2, byrow = TRUE),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20")))
    )
    d <- nmr_data(ds)
    expect_equal(rownames(d), c("10", "20"))
    expect_equal(colnames(d), c("1", "2", "3"))
    expect_equal(unname(d), matrix(c(1, 2, 1, 4, 5, 6), nrow = 2, byrow = TRUE))
})

test_that("nmr_data.nmr_dataset_1D reads a different data_ field via the what argument", {
    ds <- new_nmr_dataset_1D(
        ppm_axis = c(1, 2, 3),
        data_1r = matrix(c(1, 2, 1), nrow = 1),
        metadata = list(external = data.frame(NMRExperiment = "10"))
    )
    ds$data_1i <- matrix(c(5, 6, 7), nrow = 1)
    d <- nmr_data(ds, what = "data_1i")
    expect_equal(unname(d), matrix(c(5, 6, 7), nrow = 1))
})

test_that("nmr_data.nmr_dataset_peak_table names the rows from the sample names", {
    peak_table <- matrix(1:6, nrow = 2, dimnames = list(NULL, c("P1", "P2", "P3")))
    metadata <- list(external = data.frame(NMRExperiment = letters[1:2]))
    dataset_peak_table <- new_nmr_dataset_peak_table(peak_table, metadata)

    pt <- nmr_data(dataset_peak_table)

    expect_equal(rownames(pt), c("a", "b"))
    expect_equal(colnames(pt), c("P1", "P2", "P3"))
    expect_equal(unname(pt), unname(peak_table))
})

test_that("nmr_data has no method for unsupported classes", {
    expect_error(nmr_data(list()), "no applicable method")
})

test_that("nmr_data<-.nmr_dataset_1D replaces data_1r and strips any row/col names", {
    ds <- new_nmr_dataset_1D(
        ppm_axis = c(1, 2, 3),
        data_1r = matrix(c(1, 2, 1), nrow = 1),
        metadata = list(external = data.frame(NMRExperiment = "10"))
    )
    replacement <- matrix(c(9, 8, 7), nrow = 1, dimnames = list("10", c("a", "b", "c")))

    nmr_data(ds) <- replacement

    expect_equal(unname(ds$data_1r), matrix(c(9, 8, 7), nrow = 1))
    expect_null(rownames(ds$data_1r))
    expect_null(colnames(ds$data_1r))
})

test_that("nmr_data<-.nmr_dataset_1D requires the replacement to have one row per sample", {
    ds <- new_nmr_dataset_1D(
        ppm_axis = c(1, 2, 3),
        data_1r = matrix(c(1, 2, 1), nrow = 1),
        metadata = list(external = data.frame(NMRExperiment = "10"))
    )
    expect_error(
        nmr_data(ds) <- matrix(1:6, nrow = 2),
        "nrow"
    )
})

test_that("nmr_data<-.nmr_dataset_1D removes the field entirely when value is NULL", {
    ds <- new_nmr_dataset_1D(
        ppm_axis = c(1, 2, 3),
        data_1r = matrix(c(1, 2, 1), nrow = 1),
        metadata = list(external = data.frame(NMRExperiment = "10"))
    )
    nmr_data(ds) <- NULL
    expect_false("data_1r" %in% names(unclass(ds)))
})

test_that("nmr_data<- has no method for unsupported classes", {
    expect_error(`nmr_data<-`(list(), value = 1), "no applicable method")
})
