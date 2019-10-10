context("test-interpolation")

test_that("nmr_ppm_resolution works", {
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
  ppm_res <- nmr_ppm_resolution(dataset)[[1]]
  expect_true(is.numeric(ppm_res))
})


test_that("nmr_ppm_resolution works", {
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
  ppm_res <- nmr_ppm_resolution(dataset)[[1]]
  expect_true(is.numeric(ppm_res))
})


test_that("nmr_interpolate_1D works before interpolation", {
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
  expect_true(is.list(dataset[["axis"]]))
})


test_that("nmr_interpolate_1D works", {
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
  dataset <- nmr_interpolate_1D(dataset, axis = c(min = -0.1, max = 1, by = 0.02))
  expect_true(is.numeric(dataset[["axis"]]))
})


