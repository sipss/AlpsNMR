context("test-read_bruker")

test_that("nmr_read_samples_dir works", {
  dir_to_demo_dataset <- system.file("dataset-demo", package = "NIHSnmr")
  dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
  expect_equal(dataset$num_samples, 4)
})
