context("test-baseline_removal")

test_that("nmr_baseline_removal works", {
  dataset <- new_nmr_dataset_1D(ppm_axis = c(1,2,3),
                                data_1r = matrix(c(1,2,1), nrow = 1),
                                metadata = list(external = data.frame(NMRExperiment = "10")))
  dataset <- nmr_baseline_removal(dataset, lambda = 4, p = 0.02)
  expect_true(is.matrix(dataset[["data_1r"]]))
})
