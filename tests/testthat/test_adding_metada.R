context("test-adding_metadata")

test_that("nmr_meta_get works", {
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
  meta <- nmr_meta_get(dataset, groups = "external")
  expect_equal(meta[[1,1]], "10")
})

test_that("nmr_meta_get works", {
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
  MeOH_plasma_extraction_xlsx <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
  exp_subj_id <- readxl::read_excel(MeOH_plasma_extraction_xlsx, sheet = 1)
  dataset <- nmr_meta_add(dataset, metadata = exp_subj_id, by = "NMRExperiment")
  meta <- nmr_meta_get(dataset, groups = "external")
  expect_match(meta[[1,2]], "Ana", ignore.case = TRUE)
})

