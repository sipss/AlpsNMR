test_that("new_nmr_dataset_1D", {
  ds <- new_nmr_dataset_1D(ppm_axis = c(1,2,3),
                           data_1r = matrix(c(1,2,1), nrow = 1),
                           metadata = list(external = data.frame(NMRExperiment = "10")))
  
  expect_equal(ds[["num_samples"]], 1)
})

test_that("nmr_read_samples_dir", {
  expect_error( nmr_read_samples_dir ("XXX"))
})
