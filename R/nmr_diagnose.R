#' Deprecated function
#' @param samples An nmr_dataset
#' @seealso nmr_normalize
#' @export
#' @return normalized nmr
nmr_diagnose <- function(samples) {
    .Deprecated("nmr_normalize_extra_info")
    nmr_normalize_extra_info(samples)
}
