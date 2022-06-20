test_that("nmr_read_samples_dir works", {
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
  expect_equal(dataset$num_samples, 3)
})

test_that("nmr_read_samples returns unique NMR experiments", {
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  dataset <- nmr_read_samples(
    c(file.path(dir_to_demo_dataset, "10.zip"),
      file.path(dir_to_demo_dataset, "10.zip"))
  )
  expect_equal(dataset$num_samples, 2)
  expect_false(any(duplicated(names(dataset))))
})

test_that("create_sample_names returns good unique guesses", {
    sample_names <- c("a", "b")
    expect_equal(create_sample_names(sample_names), sample_names)
    sample_names <- c("a.zip", "b.zip")
    expect_equal(create_sample_names(sample_names), c("a", "b"))
    sample_names <- c("bar/a.zip", "bar/b.zip")
    expect_equal(create_sample_names(sample_names), c("a", "b"))
    sample_names <- c("bar/a.zip", "foo/b.zip")
    expect_equal(create_sample_names(sample_names), c("a", "b"))
    sample_names <- c("bar/a.zip", "foo/a.zip")
    expect_equal(create_sample_names(sample_names), c("bar/a", "foo/a"))
})