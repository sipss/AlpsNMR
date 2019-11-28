#' Deprecated function
#' @seealso nmr_normalize
#' @export
nmr_diagnose <- function(samples) {
  .Deprecated("nmr_normalize_extra_info")
  nmr_normalize_extra_info(samples)
}
