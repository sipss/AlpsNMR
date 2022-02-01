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
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' validate_nmr_dataset_family(dataset_1D)
validate_nmr_dataset_family <- function(nmr_dataset_family) {
    abort_if_not(
        inherits(nmr_dataset_family, "nmr_dataset_family"),
        message = "Not an nmr_dataset_family object"
    )
    abort_if_not(
        is.list(nmr_dataset_family),
        message = "nmr_dataset_family objects are list-like. This object is not"
    )
    
    num_samples <- nmr_dataset_family[["num_samples"]]
    
    abort_if_not(
        "metadata" %in% names(nmr_dataset_family),
        message = "Missing acquisition and parameter metadata"
    )

    metadata <- nmr_dataset_family[["metadata"]]
    abort_if_not(
        is.vector(metadata) & is.list(metadata),
        message = "metadata should be a list"
    )
    abort_if_not(
        "external" %in% names(metadata),
        message = "$metadata$external should be a data frame"
    )
    abort_if_not(
        all(purrr::map_lgl(metadata, is.data.frame)),
        message = "all metadata elements should be data frames"
    )
    for (metad_idx in seq_along(metadata)) {
        metad_name <- names(metadata)[metad_idx]
        metad <- metadata[[metad_idx]]
        abort_if_not(
            nrow(metad) == num_samples,
            message = glue::glue(
                "The number of rows of {metad_name} does not match the number of samples"
            )
        )
        abort_if_not(
            "NMRExperiment" %in% colnames(metad),
            message = glue::glue_data(
                list(metad_name = metad_name),
                "metadata '{metad_name}' does not include the NMRExperiment column"
            )
        )
        abort_if_not(
            all(metad[["NMRExperiment"]] == metadata[[1]][["NMRExperiment"]]),
            message = glue::glue(
                "The NMRExperiment column in {metad_name} is not equal the same column in {names(metadata)[1]}"
            )
        )
    }
    
    nmr_dataset_family
}

#' Keep samples based on metadata column criteria
#'
#' @param .data An [nmr_dataset_family] object
#' @param ... conditions, as in [dplyr]
#' @return The same object, with the matching rows
#' @importFrom dplyr filter
#' @family nmr_dataset_family manipulation functions
#' @family subsetting functions
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' 
#' ## example 1
#' sample_10 <- filter(dataset_1D, NMRExperiment == "10")
#'
#' ## example 2
#' #test_samples <- dataset_1D %>% filter(nmr_peak_table$metadata$external$Group == "placebo")
#' @export
filter.nmr_dataset_family <- function(.data, ...) {
    dots <- rlang::quos(...)
    meta <- nmr_meta_get(.data)
    meta$tmp_row_idx <- seq_len(nrow(meta))
    indices_to_keep <- dplyr::filter(meta,!!!dots)$tmp_row_idx
    return(.data[indices_to_keep])
}
