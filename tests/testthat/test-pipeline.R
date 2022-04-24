test_that("pipe_load_samples  works", {
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  expect_error(dataset <- pipe_load_samples (dir_to_demo_dataset))
})
