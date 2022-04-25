test_that("nmr_meta_get works", {
  dataset <- new_nmr_dataset_1D(ppm_axis = c(0:10),
                                data_1r = matrix(sample(0:43,replace = TRUE), nrow = 4),
                                metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30","40"))))
  dataset[["metadata"]][["external"]][["NMRExperiment"]] = as.character(dataset[["metadata"]][["external"]][["NMRExperiment"]])
  meta <- nmr_meta_get(dataset, groups = "external")
  expect_equal(meta[[1,1]], "10")
})

test_that("nmr_meta_add_tidy_excel works", {
  dataset <- new_nmr_dataset_1D(ppm_axis = c(0:10),
                                data_1r = matrix(sample(0:43,replace = TRUE), nrow = 4),
                                metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30","40"))))
  dataset[["metadata"]][["external"]][["NMRExperiment"]] = as.character(dataset[["metadata"]][["external"]][["NMRExperiment"]])
  
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  MeOH_plasma_extraction_xlsx <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
  
  dataset <- nmr_meta_add_tidy_excel(dataset, MeOH_plasma_extraction_xlsx)
  expect_match(dataset[["metadata"]][["external"]][["SubjectID"]][[1]], "Ana", ignore.case = TRUE)
})


test_that("nmr_meta_add works", {
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
  MeOH_plasma_extraction_xlsx <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
  exp_subj_id <- readxl::read_excel(MeOH_plasma_extraction_xlsx, sheet = 1)
  dataset <- nmr_meta_add(dataset, metadata = exp_subj_id, by = "NMRExperiment")
  meta <- nmr_meta_get(dataset, groups = "external")
  expect_match(meta[[1,2]], "Ana", ignore.case = TRUE)
})

