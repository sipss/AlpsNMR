## Prepare demo dataset
prepare_dataset <- function() {
  # 12 artificial samples created based on the 3 demo samples
  MeOH_plasma_extraction_dir <- system.file("dataset-demo", package = "AlpsNMR")
  MeOH_plasma_extraction_xlsx <- file.path(MeOH_plasma_extraction_dir, "dummy_metadata.xlsx")
  exp_subj_id <- readxl::read_excel(MeOH_plasma_extraction_xlsx, sheet = 1)

  zip_files <- fs::dir_ls(MeOH_plasma_extraction_dir, glob = "*.zip")
  
  dataset <- nmr_read_samples(sample_names = zip_files)
  dataset <- nmr_meta_add(dataset, metadata = exp_subj_id, by = "NMRExperiment")
  dataset <- nmr_interpolate_1D(dataset, axis = c(min = 3.7, max = 4.5, by = 2.3E-4))
  dataset <- nmr_baseline_removal(dataset, lambda = 6, p = 0.01)
  dataset <- nmr_normalize(dataset, method = "area")

  metadata <- nmr_meta_get(dataset, groups = "external")
  metadata$Group <- c("A", "B", "B")
  # Artificially create a larger dataset
  larger_metadata <- rbind(metadata, metadata, metadata, metadata, metadata)
  
  larger_metadata$NMRExperiment <- as.character(
    seq(from = 10, by = 10, length.out = nrow(larger_metadata))
  )
  data_matrix <- nmr_data(dataset)
  dataset <- new_nmr_dataset_1D(
    ppm_axis = dataset$axis,
    data_1r = rbind(data_matrix, data_matrix, data_matrix, data_matrix, data_matrix),
    metadata = list(external = larger_metadata)
  )
  dataset
}

##

is_mixomics_broken <- function() {
  # Issue https://github.com/mixOmicsTeam/mixOmics/pull/199 makes mixOmics-dependant test to fail.
  # If the issue is fixed and released, depend on that mixomics version and remove this function
  r_version <- as.package_version(paste0(R.Version()[c("major","minor")], collapse = "."))
  if (r_version < "4.2") {
    return(FALSE)
  }
  
  # test for this mixOmics issue:
  dummy_object <- list()
  class(dummy_object) <- c("mint.block.pls", "mixo_pls")
  errmsg <- tryCatch({
      mixOmics::plotIndiv(dummy_object)
    }, error = function(e)  {
      conditionMessage(e)
    }
  )
  grepl("length(class2) == 1L is not TRUE", errmsg)
}


## Dataset can be used

test_that("nmr_data_analysis works", {
  skip_if(is_mixomics_broken(), message = "Skipping nmr_data_analysis because mixOmics is broken. See https://github.com/mixOmicsTeam/mixOmics/pull/199")
  dataset <- prepare_dataset()
  methodology <- plsda_auroc_vip_method(ncomp = 2)
  set.seed(123L)
  out <- nmr_data_analysis(
    dataset,
    y_column = "Group",
    identity_column = NULL,
    external_val = list(iterations = 1, test_size = 0.25),
    internal_val = list(iterations = 2, test_size = 0.25),
    data_analysis_method = methodology
  )
  expect_false(is.null(out))
})

test_that("random subsampling works", {
  subject_id <- rep(c("Alice", "Bob", "Charlie", "Diana"), times = 2)
  replicate <- rep(c(1,2), each = 4)
  set.seed(2563432L)
  sample_idx <- 1:8
  num_iterations <- 2L
  out <- random_subsampling(sample_idx, iterations = num_iterations, test_size = 0.25,
                            keep_together = subject_id)
  expect_equal(length(out), num_iterations)
  expect_equal(length(out[[1]][["training"]]), 6L)
  expect_equal(length(out[[1]][["test"]]), 2L)
  # Subjects kept together in the split, no subject in train is present in test:
  expect_equal(
    length(
      intersect(
        subject_id[out[[1]][["test"]]],
        subject_id[out[[1]][["training"]]]
      )
    ),
    0L
  )
})

test_that("split_double_cv works", {
  nsamples <- 16L
  subject_id <- rep(c("Alice", "Bob", "Charlie", "Diana"), times = 4)
  replicate <- rep(c(1,2), each = 8)
  metadata <- data.frame(
    NMRExperiment = as.character(seq(from = 10, by = 10, length.out = nsamples)),
    SubjectID = subject_id,
    Replicate = replicate,
    stringsAsFactors = FALSE
  )
  dataset <- new_nmr_dataset_1D(
    ppm_axis = 1:10,
    data_1r = matrix(sample(1:200, 10*nsamples), ncol = 10, nrow = nsamples),
    metadata = list(external = metadata)
  )
  
  external_val_niter <- 2L
  internal_val_niter <- 4L
  external_test_size <- 0.25
  internal_test_size <- 0.34
  out <- AlpsNMR:::split_double_cv(
    dataset = dataset,
    keep_together = "SubjectID",
    external_val = list(iterations = external_val_niter, test_size = external_test_size),
    internal_val = list(iterations = internal_val_niter, test_size = internal_test_size)
  )
  
  expect_equal(names(out), c("outer", "inner"))
  expect_equal(length(out[["outer"]]), external_val_niter)
  expect_equal(length(out[["inner"]]), external_val_niter*internal_val_niter)
  expected_samples_in_external_test <- floor(nsamples*external_test_size)
  expected_samples_in_train <- nsamples - expected_samples_in_external_test
  expected_samples_in_train_internal_test <- floor(expected_samples_in_train*internal_test_size)
  expected_samples_in_train_internal_train <- expected_samples_in_train - expected_samples_in_train_internal_test
  
  expect_equal(length(out$inner$`1_1`$inner_train_idx),
               expected_samples_in_train_internal_train)
})
