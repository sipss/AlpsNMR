#' Validate 1D nmr datasets
#' @param nmr_dataset_1D An nmr_dataset_1D object
#' @return The nmr_dataset_1D unchanged
#' 
#' This function is useful for its side-effects. Stopping in case of error
#' 
#' @export
validate_nmr_dataset_1D <- function(nmr_dataset_1D) {
  assert_that(inherits(nmr_dataset_1D, "nmr_dataset_1D"),
              msg = "Not an nmr_dataset_1D")
  assert_that(is.list(nmr_dataset_1D),
              msg = "nmr_dataset_1D objects are list-like. This object is not")
  
  assert_that("axis" %in% names(nmr_dataset_1D),
              msg = "nmr_dataset_1D must have a ppm axis")
  assert_that("data_1r" %in% names(nmr_dataset_1D),
              msg = "nmr_dataset_1D must have a data_1r matrix")
  
  assert_that("metadata" %in% names(nmr_dataset_1D), msg = "Missing acquisition and parameter metadata")
  assert_that("metadata_ext" %in% names(nmr_dataset_1D), msg = "Missing external metadata")
  metadata <- nmr_dataset_1D[["metadata"]]
  metadata_ext <- nmr_dataset_1D[["metadata_ext"]]
  assert_that(is.data.frame(metadata), msg = "metadata should be a data frame")
  assert_that(is.data.frame(metadata_ext), msg = "metadata_ext should be a data frame")
  assert_that("NMRExperiment" %in% colnames(metadata), msg = "metadata should have an NMRExperiment column")
  assert_that("NMRExperiment" %in% colnames(metadata_ext), msg = "metadata_ext should have an NMRExperiment column")
  assert_that(nrow(metadata) == nrow(metadata_ext), msg = "metadata and metadata_ext should have the same number of rows")
  
  ppm_axis <- nmr_dataset_1D[["axis"]]
  data_1r <- nmr_dataset_1D[["data_1r"]]
  assert_that(is.vector(ppm_axis) && is.numeric(ppm_axis),
              msg = "axis must be a numeric vector")
  assert_that(is.matrix(data_1r) && is.numeric(data_1r),
              msg = "data_1r must be a numeric matrix")
  
  assert_that(length(ppm_axis) == ncol(data_1r),
              msg = "ppm axis does not have a length equal to ncol(data_1r)")
  num_samples <- nrow(data_1r)
  
  assert_that(num_samples == nmr_dataset_1D[["num_samples"]],
              msg = "The num_samples value does not match nrow(data_1r)")
  assert_that(num_samples == nrow(metadata),
              msg = "data_1r and acquisition metadata differ in number of samples")
  assert_that(nrow(metadata_ext) == nrow(metadata),
              msg = "Acquisition metadata and external metadata differ in number of rows")
  nmr_dataset_1D
}

#' Creates a new 1D nmr_dataset object from scratch
#' 
#' @param ppm_axis A numeric vector with the ppm values for the columns of data_1r
#' @param data_1r A numeric matrix with one NMR spectrum on each row
#' @param metadata A data frame with at least the `NMRExperiment` column
#' @param metadata_ext A data frame with at least the `NMRExperiment` column
#' 
#' @importFrom assertthat assert_that
#' @importFrom glue glue
#' @export
new_nmr_dataset_1D <- function(ppm_axis, data_1r, metadata, metadata_ext) {
  samples <- list()
  samples[["metadata_ext"]] <- metadata_ext
  samples[["metadata"]] <- metadata
  samples[["data_1r"]] <- data_1r
  samples[["axis"]] <- ppm_axis
  samples[["num_samples"]] <- nrow(data_1r)
  samples[["processing"]] <- list(data_loaded = TRUE,
                                  interpolation = TRUE,
                                  exclusion = FALSE,
                                  normalization = FALSE)
  class(samples) <- "nmr_dataset_1D"
  validate_nmr_dataset_1D(samples)
  samples
}
