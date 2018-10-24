#' Exclude region from samples
#'
#' Sets to zero a given region (for instance to remove the water peak)
#' @param samples A [nmr_dataset] object
#' @param exclude A list with regions to be zeroed. Typically:
#'                `exclude = list(water = c(4.7, 5.0))`
#' @return The [nmr_dataset] object, with the regions excluded
#' @export
nmr_exclude_region <- function(samples, exclude = list(water = c(4.7, 5.0))) {
  if (is.null(exclude) || length(exclude) == 0) {
    return(samples)
  }
  if (!samples[["processing"]][["interpolation"]]) {
    stop("exclusion not implemented for non-interpolated samples")
  }
  
  samples[["processing"]][["exclusion_params"]] <- exclude
  data_fields <- names(samples)[grepl(pattern = "^data_.*", x = names(samples))]
  for (data_field in data_fields) {
    dimension <- length(dim(samples[[data_field]])) - 1 # first dim is samples
    for (i_region in seq_along(exclude)) {
      region <- exclude[[i_region]]
      if (dimension == 1) {
        
        excl_dim1 <- samples[["axis"]][[1]] >= min(region) & samples[["axis"]][[1]] <= max(region)
        samples[[data_field]][,excl_dim1] <- 0
      } else if (dimension == 2) {
        excl_dim1 <- samples[["axis"]][[1]] >= min(region[1:2]) & samples[["axis"]][[1]] <= max(region[1:2])
        excl_dim2 <- samples[["axis"]][[2]] >= max(region[3:4]) & samples[["axis"]][[2]] <= max(region[3:4])
        samples[[data_field]][,excl_dim1, excl_dim2] <- 0
      } else {
        stop("Not implemented error")
      }
    }
  }
  samples[["processing"]][["exclusion"]] <- TRUE
  return(samples)
}