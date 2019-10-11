context("test-integration")

test_that("nmr_integrate_regions works", {
  dataset <- new_nmr_dataset_1D(ppm_axis = 1:10,
                                data_1r = matrix(sample(0:99,replace = TRUE), nrow = 10),
                                metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100"))))
  peak_table_integration = nmr_integrate_regions(
  samples = dataset,
  regions = list(random = c(2,5)),
  fix_baseline = TRUE)
  expect_true(is.matrix(peak_table_integration[["peak_table"]]))
  expect_true(is.numeric(peak_table_integration[["peak_table"]][[1]]))
})

test_that("nmr_integrate_peak_positions works", {
  dataset <- new_nmr_dataset_1D(ppm_axis = 1:10,
                                data_1r = matrix(sample(0:99,replace = TRUE), nrow = 10),
                                metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100"))))
  peak_table_integration = nmr_integrate_peak_positions(
  samples = dataset,
  peak_pos_ppm = list(c(2,3,4)),
  peak_width_ppm = 1)
  expect_true(is.matrix(peak_table_integration[["peak_table"]]))
  expect_true(is.numeric(peak_table_integration[["peak_table"]][[1]]))
})


