test_that("nmr_exclude_region works", {
  dataset <- new_nmr_dataset_1D(ppm_axis = c(1,2,3),
                           data_1r = matrix(c(1,2,1), nrow = 1),
                           metadata = list(external = data.frame(NMRExperiment = "10")))
  regions_to_exclude <- list(water = c(1, 2))
  dataset <- nmr_exclude_region(dataset, exclude = regions_to_exclude)
  expect_equal(dataset[["data_1r"]], matrix(1, nrow = 1))
})

test_that("AlpsNMR filter works", {
  dataset <- new_nmr_dataset_1D(ppm_axis = c(1,2,3),
                                data_1r = matrix(c(1,2,1,1,2,1), nrow = 2),
                                metadata = list(external = data.frame(NMRExperiment = c("10", "20"))))
  dataset <- filter(dataset, NMRExperiment == "10")
  expect_equal(dataset[["num_samples"]], 1)
})

test_that("AlpsNMR filter works", {
  dataset <- new_nmr_dataset_1D(ppm_axis = c(1,2,3),
                                data_1r = matrix(c(1,2,1,1,2,1), nrow = 2),
                                metadata = list(external = data.frame(NMRExperiment = c("10", "20"))))
  dataset <- filter(dataset, NMRExperiment == "10")
  meta <- nmr_meta_get(dataset, groups = "external")
  #Force error language to English
  Sys.setenv(LANGUAGE='en')
  expect_error(meta[["NMRExperiment"]][[2]], "subscript out of bounds")
})
