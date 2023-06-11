test_that("alps_asics works", {
    skip_if_not_installed("ASICS")
    dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
    dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
    dataset <- nmr_interpolate_1D(dataset, axis = c(min = 1, max = 2, by = 0.002))
    spec_obj <- alps_asics(dataset)
    expect_true(class(spec_obj)[1]=="Spectra")
})