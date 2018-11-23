# Needed for filter.nmr_dataset:
#' @importFrom dplyr filter
#' @export
dplyr::filter

#' Keep samples based on metadata column criteria
#'
#' @param .data An [nmr_dataset] object
#' @param ... conditions, as in [dplyr]
#' @return The same object, with the matching rows
#' @importFrom dplyr filter
#' @family nmr_dataset manipulation functions
#' @family subsetting functions
#' @export
filter.nmr_dataset <- function(.data, ...) {
  dots <- rlang::quos(...)
  meta <- nmr_meta_get(.data)
  meta$tmp_row_idx <- seq_len(nrow(meta))
  indices_to_keep <- dplyr::filter(meta, !!! dots)$tmp_row_idx
  return(.data[indices_to_keep])
}

#' @rdname filter.nmr_dataset
#' @family nmr_dataset_1D manipulation functions
#' @family subsetting functions
#' @export
filter.nmr_dataset_1D <- filter.nmr_dataset

#' @rdname filter.nmr_dataset
#' @family nmr_dataset_peak_table manipulation functions
#' @family subsetting functions
#' @export
filter.nmr_dataset_peak_table <- filter.nmr_dataset
