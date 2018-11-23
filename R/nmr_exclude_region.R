#' Exclude region from samples
#'
#' Excludes a given region (for instance to remove the water peak)
#' 
#' @param samples An object
#' @param exclude A list with regions to be removed Typically:
#'                `exclude = list(water = c(4.7, 5.0))`
#' @return The same object, with the regions excluded
#' @export
nmr_exclude_region <- function(samples, exclude = list(water = c(4.7, 5.0))) {
  UseMethod("nmr_exclude_region")
}

#' @rdname nmr_exclude_region
#' @family nmr_dataset_1D functions
#' @export
nmr_exclude_region.nmr_dataset_1D <- function(samples, exclude = list(water = c(4.7, 5.0))) {
  if (is.null(exclude) || length(exclude) == 0) {
    return(samples)
  }
  axis_include <- rep(TRUE, length(samples[["axis"]]))
  for (i_region in seq_along(exclude)) {
    region <- exclude[[i_region]]
    excl_dim1 <- samples[["axis"]] >= min(region) & samples[["axis"]] <= max(region)
    axis_include[excl_dim1] <- FALSE
  }
  samples[["axis"]] <- samples[["axis"]][axis_include]
  samples[["data_1r"]] <- samples[["data_1r"]][, axis_include]
  return(samples)
}

