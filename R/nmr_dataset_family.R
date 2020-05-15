#' nmr_dataset like objects (S3 classes)
#'
#' The AlpsNMR package defines and uses several objects to manage NMR Data.
#'
#' These objects share some structure and functions, so it makes sense to have
#' an abstract class to ensure that the shared structures are compatible
#'
#' @name nmr_dataset_family
#' @family AlpsNMR dataset objects
#' @seealso [Functions to save and load these objects][load_and_save_functions]
NULL



#' Validate nmr_dataset_family objects
#' @param nmr_dataset_family An [nmr_dataset_family] object
#' @return The [nmr_dataset_family] unchanged
#'
#' This function is useful for its side-effects: Stopping in case of error
#'
#' @family nmr_dataset_family functions
#' @family class helper functions
#' @export
validate_nmr_dataset_family <- function(nmr_dataset_family) {
    assert_that(inherits(nmr_dataset_family, "nmr_dataset_family"),
                            msg = "Not an nmr_dataset_family object")
    assert_that(is.list(nmr_dataset_family),
                            msg = "nmr_dataset_family objects are list-like. This object is not")
    
    num_samples <- nmr_dataset_family[["num_samples"]]
    
    assert_that("metadata" %in% names(nmr_dataset_family), msg = "Missing acquisition and parameter metadata")
    metadata <- nmr_dataset_family[["metadata"]]
    assert_that(is.vector(metadata) &
                                is.list(metadata), msg = "metadata should be a list")
    assert_that("external" %in% names(metadata),
                            msg = "$metadata$external should be a data frame")
    assert_that(all(purrr::map_lgl(metadata, is.data.frame)), msg = "all metadata elements should be data frames")
    for (metad_idx in seq_along(metadata)) {
        metad_name <- names(metadata)[metad_idx]
        metad <- metadata[[metad_idx]]
        assert_that(
            nrow(metad) == num_samples,
            msg = glue::glue(
                "The number of rows of {metad_name} does not match the number of samples"
            )
        )
        assert_that(
            "NMRExperiment" %in% colnames(metad),
            msg = glue::glue_data(
                list(metad_name = metad_name),
                "metadata '{metad_name}' does not include the NMRExperiment column"
            )
        )
        assert_that(
            all(metad[["NMRExperiment"]] == metadata[[1]][["NMRExperiment"]]),
            msg = glue::glue(
                "The NMRExperiment column in {metad_name} is not equal the same column in {names(metadata)[1]}"
            )
        )
    }
    
    nmr_dataset_family
}



#' @importFrom dplyr filter
#' @export
dplyr::filter

#' Keep samples based on metadata column criteria
#'
#' @param .data An [nmr_dataset_family] object
#' @param ... conditions, as in [dplyr]
#' @return The same object, with the matching rows
#' @importFrom dplyr filter
#' @family nmr_dataset_family manipulation functions
#' @family subsetting functions
#' @examples
#' \dontrun{
#' ## example 1
#' placebo_samples <- filter(nmr_dataset, Group == "placebo")
#'
#' ## example 2
#' test_samples <- nmr_dataset %>% filter(nmr_peak_table$metadata$external$Group == "placebo")
#' }
#' @export
filter.nmr_dataset_family <- function(.data, ...) {
    dots <- rlang::quos(...)
    meta <- nmr_meta_get(.data)
    meta$tmp_row_idx <- seq_len(nrow(meta))
    indices_to_keep <- dplyr::filter(meta,!!!dots)$tmp_row_idx
    return(.data[indices_to_keep])
}
