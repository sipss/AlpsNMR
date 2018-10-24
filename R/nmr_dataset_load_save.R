#' Load a [nmr_dataset] from a file
#'
#' Loads a [nmr_dataset] saved with [nmr_dataset_save()]
#'
#' @param file_name The file name with the [nmr_dataset]
#' @return the [nmr_dataset] object stored in `file_name`
#' @export
nmr_dataset_load <- function(file_name) {
  return(readRDS(file_name))
}

#' Save an [nmr_dataset] object
#'
#' Wraps `save` to save the [nmr_dataset] object in `.RDS` format
#'
#' @param nmr_dataset The [nmr_dataset] object to save
#' @param file_name The output file name.
#' @param ... Optional arguments passed to [saveRDS].
#' @return the passed [nmr_dataset] object
#' @export
nmr_dataset_save <- function(nmr_dataset, file_name, ...) {
  saveRDS(nmr_dataset, file_name)
  return(nmr_dataset)
}

# This function can be removed once we know it is not needed by our users
nmr_dataset_load_old_and_save <- function(old_file_name, new_file_name) {
  nmr_dataset <- NULL
  load(old_file_name)
  if (is.null(nmr_dataset)) {
    stop("This was not saved with the NIHSnmr version 1.0 or lower")
  }
  nmr_dataset_save(nmr_dataset, new_file_name)
}