test_that("nmr_ppm_resolution works", {
  dataset <- new_nmr_dataset_1D(ppm_axis = c(1,2,3),
                                data_1r = matrix(c(1,2,1), nrow = 1),
                                metadata = list(external = data.frame(NMRExperiment = "10")))
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


