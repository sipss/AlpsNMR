test_that("nmr_pca_outliers_robust works", {
  dataset <- new_nmr_dataset_1D(ppm_axis = c(0:10),
                                data_1r = matrix(sample(0:43,replace = FALSE), nrow = 4),
                                metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30","40"))))
  dataset[["metadata"]][["external"]][["NMRExperiment"]] = as.character(dataset[["metadata"]][["external"]][["NMRExperiment"]])
  pca_outliers <- nmr_pca_outliers_robust(dataset)
  
  pca_built = nmr_pca_build_model (dataset)
  plot_variance = nmr_pca_plot_variance (pca_built)
  score = nmr_pca_scoreplot (dataset, pca_built)
  loadings = nmr_pca_loadingplot (pca_built,2)
  
  pca_outliers_no_robust = nmr_pca_outliers(dataset, pca_built)
  outliers_plot = nmr_pca_outliers_plot (dataset, pca_outliers_no_robust)
  
  plot = nmr_pca_outliers_plot(dataset, pca_outliers)
  dataset <- nmr_pca_outliers_filter(dataset, pca_outliers)
  
  
  expect_true(is.numeric(pca_outliers[["outlier_info"]][["Tscores"]]))
  expect_true(is.numeric(pca_outliers_no_robust[["outlier_info"]][["Tscores"]]))
  
  expect_true(is.matrix(pca_built[["X"]]))
  expect_true(is.list(plot_variance))
  expect_true(is.list(score))
  expect_true(is.list(loadings))
  expect_true(is.list(outliers_plot))
  expect_true(is.list(plot))
  expect_true(is.integer(dataset[["data_1r"]][[1]]))
})
