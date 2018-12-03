#' Set/Get diagnostic information
#' 
#' Set or get diagnostic information
#'
#' @param nmr_dataset An [nmr_dataset_family] object
#'
#' @return Diagnostic information, usually a list with data or plots
#' @export
#'
nmr_diagnose <- function(nmr_dataset) {
  attr(nmr_dataset, "diagnostic")
}


#' @rdname nmr_diagnose
#' @param value The diagnostic we want to set
#' @export
"nmr_diagnose<-" <- function(nmr_dataset, value) {
  attr(nmr_dataset, "diagnostic") <- value
  nmr_dataset
}