test_that("to_ChemoSpec works", {
  skip_if_not_installed("ChemoSpec")
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
  dataset <- nmr_interpolate_1D(dataset, axis = c(min = 1, max = 2, by = 0.002))
  Ch <- to_ChemoSpec(dataset)
  expect_true(is.list(Ch))
})

test_that("to_ChemoSpec works with a group", {
  skip_if_not_installed("ChemoSpec")
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
  dataset <- nmr_interpolate_1D(dataset, axis = c(min = 1, max = 2, by = 0.002))
  extra_metadata <- data.frame(
    NMRExperiment = c("10", "20", "30"),
    treatment = c("A", "B", "A")
  )
  dataset <- nmr_meta_add(dataset, extra_metadata)
  Ch <- to_ChemoSpec(dataset, group = "treatment")
  expect_true(is.list(Ch))
  xx <- ChemoSpec::sumGroups(Ch)
  expect_equal(xx$group, c("A", "B"))
  expect_equal(xx$no., c(2L, 1L))
})

