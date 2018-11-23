#' nmr_dataset_peak_table (S3 class)
#'
#' An `nmr_dataset_peak_table` represents a peak table with metadata.
#' It is defined as an S3 class, and it can be treated as a regular list.
#'
#' - `metadata`: A list of data frames. Each data frame contains metadata. Usually
#' the list only has one data frame named "external".
#' 
#' - `peak_table`: A matrix with one sample on each row and the peaks in the
#' columns
#' 
#' @name nmr_dataset_peak_table
NULL

#' Validate nmr_dataset_peak_table objects
#' @param nmr_dataset_peak_table An [nmr_dataset_peak_table] object
#' @return The [nmr_dataset_peak_table] unchanged
#' 
#' This function is useful for its side-effects: Stopping in case of error
#' 
#' @export
validate_nmr_dataset_peak_table <- function(nmr_dataset_peak_table) {
  assert_that(inherits(nmr_dataset_peak_table, "nmr_dataset_peak_table"),
              msg = "Not an nmr_dataset_peak_table")
  assert_that(is.list(nmr_dataset_peak_table),
              msg = "nmr_dataset_peak_table objects are list-like. This object is not")
  
  assert_that("peak_table" %in% names(nmr_dataset_peak_table),
              msg = "nmr_dataset_peak_table must have a peak_table matrix")
  
  peak_table <- nmr_dataset_peak_table[["peak_table"]]
  assert_that(is.matrix(peak_table) && is.numeric(peak_table),
              msg = "peak_table must be a numeric matrix")
  
  num_samples <- nrow(peak_table)
  
  assert_that(num_samples == nmr_dataset_peak_table[["num_samples"]],
              msg = "The num_samples value does not match nrow(peak_table)")
  
  
  assert_that("metadata" %in% names(nmr_dataset_peak_table), msg = "Missing acquisition and parameter metadata")
  metadata <- nmr_dataset_peak_table[["metadata"]]
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
  nmr_dataset_peak_table
}

#' Creates a new nmr_dataset_peak_table object from scratch
#' 
#' @param peak_table A numeric matrix with one NMR spectrum on each row
#' @param metadata A list of data frames with at least the `NMRExperiment` column
#' 
#' @importFrom assertthat assert_that
#' @importFrom glue glue
#' @export
new_nmr_dataset_peak_table <- function(peak_table, metadata) {
  samples <- list()
  samples[["metadata"]] <- metadata
  samples[["peak_table"]] <- peak_table
  samples[["num_samples"]] <- nrow(peak_table)
  class(samples) <- "nmr_dataset_peak_table"
  validate_nmr_dataset_peak_table(samples)
  samples
}

#' Object is of [nmr_dataset_peak_table] class
#' @param x An object
#' @return `TRUE` if the object is an `nmr_dataset_peak_table`, `FALSE` otherwise
#' @export
is.nmr_dataset_peak_table <- function(x) inherits(x, "nmr_dataset_peak_table")

#' @export
print.nmr_dataset_peak_table <- function(x, ...) {
  cat(format(x, ...), "\n")
  invisible(x)
}

#' @export
format.nmr_dataset_peak_table <- function(x, ...) {
  paste0("An nmr_dataset_peak_table (", x$num_samples, " samples, and ", ncol(x$peak_table), " peaks)")
}

#' Extract parts of an nmr_dataset_peak_table
#' @param x an [nmr_dataset_peak_table] object
#' @param i indices of the samples to keep
#' @return an nmr_dataset_peak_table with the extracted samples
#' @export
`[.nmr_dataset_peak_table` <- function(x, i) {
  output <- x
  output$metadata <- purrr::map(output$metadata, function(metad) {
    metad[i, , drop = FALSE]
  })
  output[["peak_table"]] <- output[["peak_table"]][i, , drop = FALSE]
  output$num_samples <- nrow(output$metadata[[1]])
  validate_nmr_dataset_peak_table(output)
  return(output)
}
