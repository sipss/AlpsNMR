context("test-to_ChemoSpec")

test_that("to_ChemoSpec works", {
  skip_if_not_installed("ChemoSpec")
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
  data2 <- dataset[1:3]
  dataset <- nmr_interpolate_1D(dataset, axis = c(min = 1, max = 2, by = 0.002))
  Ch <- to_ChemoSpec(dataset)
  expect_true(is.list(Ch))
})
