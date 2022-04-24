test_that("nmr_baseline_removal works", {
  dataset <- new_nmr_dataset_1D(ppm_axis = c(1,2,3,4),
                                data_1r = matrix(c(1,2,2,1), nrow = 1),
                                metadata = list(external = data.frame(NMRExperiment = "10")))
  dataset <- nmr_baseline_removal(dataset, lambda = 4, p = 0.02)
  expect_true(is.matrix(dataset[["data_1r"]]))
})

test_that("nmr_baseline_threshold works", {
  dataset <- new_nmr_dataset_1D(ppm_axis = c(9.5, 9.7, 10.0),
                                data_1r = matrix(sample(0:2, 6, replace = TRUE), nrow = 2),
                                metadata = list(external = data.frame(NMRExperiment = c("10","20"))))
  
  n <- nmr_baseline_threshold (dataset)
  expect_true(is.numeric(n))
})
