new_nmr_dataset_2D <- function(f1_axis, f2_axis, data_2rr, metadata) {
  samples <- list()
  samples[["metadata"]] <- metadata
  samples[["data_2rr"]] <- data_2rr
  samples[["f1_axis"]] <- f1_axis
  samples[["f2_axis"]] <- f2_axis
  samples[["num_samples"]] <- dim(data_2rr)[1]
  samples[["processing"]] <- list(data_loaded = TRUE,
                                  interpolation = TRUE,
                                  exclusion = FALSE,
                                  normalization = FALSE)
  class(samples) <- "nmr_dataset_2D"
  validate_nmr_dataset_2D(samples)
  samples
}

validate_nmr_dataset_2D <- function(samples) {
  samples
}