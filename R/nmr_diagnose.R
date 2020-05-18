#' Deprecated function
#' @param samples An nmr_dataset
#' @seealso nmr_normalize
#' @export
#' @return normalized nmr
#' @examples 
#' nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' nmr_diagnose(nmr_dataset)
nmr_diagnose <- function(samples) {
    .Deprecated("nmr_normalize_extra_info")
    nmr_normalize_extra_info(samples)
}
