context("test-outliers")

test_that("nmr_pca_outliers_robust works", {
  dataset <- new_nmr_dataset_1D(ppm_axis = c(0:10),
                                data_1r = matrix(sample(0:43,replace = TRUE), nrow = 4),
                                metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30","40"))))
  # dataset[["metadata"]][["external"]][["NMRExperiment"]] = as.character(dataset[["metadata"]][["external"]][["NMRExperiment"]])
  pca_outliers <- nmr_pca_outliers_robust(dataset)
  dataset <- nmr_pca_outliers_filter(dataset, pca_outliers)
  expect_true(is.numeric(pca_outliers[["outlier_info"]][["Tscores"]]))
  expect_true(is.integer(dataset[["data_1r"]][[1]]))
})
