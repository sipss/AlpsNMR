context("test-normalization")

test_that("nmr_normalize & nmr_diagnose work", {
  dataset <- new_nmr_dataset_1D(ppm_axis = 1:10,
                                data_1r = matrix(sample(0:99,replace = TRUE), nrow = 10),
                                metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100"))))
  dataset <- nmr_normalize(dataset, method = "pqn")
  diagnostic <- nmr_diagnose(dataset)
  expect_true(is.matrix(dataset[["data_1r"]]))
  expect_true(is.numeric(dataset[["data_1r"]][[1]]))
  expect_true(is.list(diagnostic))
  expect_true(is.data.frame(diagnostic[["norm_factor"]]))
  expect_true(is.factor(diagnostic[["norm_factor"]][[1,1]]))
})


