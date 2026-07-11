test_that("nmr_ppm_resolution works", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = c(1, 2, 3),
        data_1r = matrix(c(1, 2, 1), nrow = 1),
        metadata = list(external = data.frame(NMRExperiment = "10"))
    )
    ppm_res <- nmr_ppm_resolution(dataset)[[1]]
    expect_true(is.numeric(ppm_res))
})


test_that("nmr_interpolate_1D works", {
    dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
    dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
    dataset_interpolated <- nmr_interpolate_1D(dataset, axis = c(min = -0.1, max = 1, by = 0.02))
    expect_true(is.list(dataset[["axis"]]))
    expect_true(is.numeric(dataset_interpolated[["axis"]]))
})

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
