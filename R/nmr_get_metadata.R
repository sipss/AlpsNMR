#' Get metadata
#' @param samples a [nmr_dataset] object
#' @param columns Columns to get. By default gets all the columns.
#' @param simplify Removes columns that are constant along all samples
#' @return a data frame with the injection metadata
#' @export
nmr_get_metadata <- function(samples, columns = NULL, simplify = FALSE) {
  metadata <- dplyr::left_join(samples[["metadata_ext"]],
                               samples[["metadata"]],
                               by = "NMRExperiment")
  
  # Default columns means all columns
  if (is.null(columns)) {
    columns <- colnames(metadata)
  }
  
  # NMRExperiment is always present in the output
  if (!"NMRExperiment" %in% columns) {
    columns <- c("NMRExperiment", columns)
  }
  
  # Report if user wants columns that are not present:
  if (!all(columns %in% colnames(metadata))) {
    cols_miss <- columns[!columns %in% colnames(metadata)]
    # If there are less than 10 missing columns warn all the missing names,
    # otherwise warn the first 5
    if (length(cols_miss) < 10) {
      show_cols <- length(cols_miss)
    } else {
      show_cols <- 5
    }
    warning("Missing columns: ",
            paste(utils::head(cols_miss, n = show_cols), collapse = ", "),
            " and ", length(cols_miss) - show_cols, " columns more.")
    rm(cols_miss, show_cols)
  }
  
  columns <- columns[columns %in% colnames(metadata)]
  # drop = FALSE ensures we never return a vector (always a data frame/tibble)
  metadata <- metadata[,columns, drop = FALSE]
  if (simplify) {
    metadata <- simplify_df(metadata)$diff
  }
  return(metadata)
}
