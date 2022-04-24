test_that("nmr_normalize & nmr_normalize_extra_info work", {
  dataset <- new_nmr_dataset_1D(
    ppm_axis = 1:10,
    data_1r = matrix(sample(0:99,replace = TRUE), nrow = 10),
    metadata = list(external = data.frame(
      NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100"),
      stringsAsFactors = FALSE
    ))
  )
  dataset <- nmr_normalize(dataset, method = "pqn")
  diagnostic <- nmr_normalize_extra_info(dataset)
  expect_true(is.matrix(dataset[["data_1r"]]))
  expect_true(is.numeric(dataset[["data_1r"]][[1]]))
  expect_true(is.list(diagnostic))
  expect_true(is.data.frame(diagnostic[["norm_factor"]]))
  expect_true(is.character(diagnostic[["norm_factor"]][[1,1]]))
})

test_that("nmr_normalize with negative values and below 10 samples work", {
  dataset <- new_nmr_dataset_1D(ppm_axis = 1:10,
                                data_1r = matrix(sample(-1:88,replace = FALSE), nrow = 9),
                                metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90"))))
  expect_warning(nmr_normalize(dataset, method = "pqn"))
})

test_that("nmr_normalize works with unknown method", {
  dataset <- new_nmr_dataset_1D(ppm_axis = 1:10,
                                data_1r = matrix(sample(-1:88,replace = FALSE), nrow = 9),
                                metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90"))))
  dataset_none <- nmr_normalize(dataset, method = "none")
  dataset_area <- nmr_normalize(dataset, method = "area")
  dataset_max <- nmr_normalize(dataset, method = "max")
  dataset_ <- nmr_normalize(dataset)
  
  expect_true(is.matrix(dataset_none[["data_1r"]]))
  expect_true(is.matrix(dataset_area[["data_1r"]]))
  expect_true(is.matrix(dataset_max[["data_1r"]]))
  expect_true(is.matrix(dataset_[["data_1r"]]))
})
