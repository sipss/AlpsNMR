#' nmr_dataset_1D (S3 class)
#'
#' An `nmr_dataset_1D` represents a set of 1D interpolated NMR samples.
#' It is defined as an S3 class, and it can be treated as a regular list.
#'
#' It currently has the following elements:
#'
#' 
#' - `metadata`: A list of data frames. Each data frame contains metadata of
#' a given area (acquisition parameters, preprocessing parameters, general sample information...)
#' 
#' - `axis`: A numeric vector with the chemical shift axis in ppm.
#' 
#' - `data_1r`: A matrix with one sample on each row and the chemical 
#' shifts in the columns.
#' 
#' 
#' @name nmr_dataset_1D
NULL

#' Validate 1D nmr datasets
#' @param nmr_dataset_1D An [nmr_dataset_1D] object
#' @return The [nmr_dataset_1D] unchanged
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
  
  
  assert_that("metadata" %in% names(nmr_dataset_1D), msg = "Missing acquisition and parameter metadata")
  metadata <- nmr_dataset_1D[["metadata"]]
  assert_that(is.vector(metadata) & is.list(metadata), msg = "metadata should be a list")
  assert_that("external" %in% names(metadata),
              msg = "$metadata$external should be a data frame")
  assert_that(all(purrr::map_lgl(metadata, is.data.frame)), msg = "all metadata elements should be data frames")
  for (metad_idx in seq_along(metadata)) {
    metad_name <- names(metadata)[metad_idx]
    metad <- metadata[[metad_idx]]
    assert_that(nrow(metad) == num_samples,
                msg = glue::glue("The number of rows of {metad_name} does not match the number of samples"))
    assert_that("NMRExperiment" %in% colnames(metad),
                msg = glue::glue_data(
                  list(metad_name = metad_name),
                  "metadata '{metad_name}' does not include the NMRExperiment column"))
    assert_that(all(metad[["NMRExperiment"]] == metadata[[1]][["NMRExperiment"]]),
                msg = glue::glue("The NMRExperiment column in {metad_name} is not equal the same column in {names(metadata)[1]}"))
  }
  nmr_dataset_1D
}

#' Creates a new 1D nmr_dataset object from scratch
#' 
#' @param ppm_axis A numeric vector with the ppm values for the columns of data_1r
#' @param data_1r A numeric matrix with one NMR spectrum on each row
#' @param metadata A list of data frames with at least the `NMRExperiment` column
#' 
#' @importFrom assertthat assert_that
#' @importFrom glue glue
#' @export
new_nmr_dataset_1D <- function(ppm_axis, data_1r, metadata) {
  samples <- list()
  samples[["metadata"]] <- metadata
  samples[["data_1r"]] <- data_1r
  samples[["axis"]] <- ppm_axis
  samples[["num_samples"]] <- nrow(data_1r)
  class(samples) <- "nmr_dataset_1D"
  validate_nmr_dataset_1D(samples)
  samples
}

#' Object is of [nmr_dataset_1D] class
#' @param x An object
#' @return `TRUE` if the object is an `nmr_dataset_1D`, `FALSE` otherwise
#' @export
is.nmr_dataset_1D <- function(x) inherits(x, "nmr_dataset_1D")

#' @export
print.nmr_dataset_1D <- function(x, ...) {
  cat(format(x, ...), "\n")
  invisible(x)
}

#' @export
format.nmr_dataset_1D <- function(x, ...) {
  paste0("An nmr_dataset_1D (", x$num_samples, " samples)")
}

#' Extract parts of an nmr_dataset_1D
#' @param x an [nmr_dataset_1D] object
#' @param i indices of the samples to keep
#' @return an nmr_dataset_1D with the extracted samples
#' @export
`[.nmr_dataset_1D` <- function(x, i) {
  output <- x
  output$metadata <- purrr::map(output$metadata, function(metad) {
    metad[i, , drop = FALSE]
  })
  output[["data_1r"]] <- output[["data_1r"]][i, , drop = FALSE]
  output$num_samples <- nrow(output$metadata[[1]])
  validate_nmr_dataset_1D(output)
  return(output)
}
