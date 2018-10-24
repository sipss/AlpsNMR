#' PQN normalization
#' @noRd
#' @param spectra A matrix with one spectrum on each row
#' @return A matrix with one spectrum on each row (normalized)
norm_pqn <- function(spectra) {
  num_samples <- nrow(spectra)
  if (num_samples < 10) {
    warning("The Probabalistic Quotient Normalization requires several samples ",
            "to compute the median spectra. Your number of samples is low")
  }
  # Normalize to the area
  areas <- rowSums(spectra)
  spectra2 <- spectra / areas
  if (num_samples == 1) {
    # We have warned, and here there is nothing to do anymore
    warning("PQN is absurd with a single sample. We have normalized it to the area.")
    return(list(spectra = spectra2,
                norm_factor = areas))
  }
  # Move spectra above zero:
  if (any(spectra2 < 0)) {
    spectra2 <- spectra2 - min(spectra2)
  }
  # Median of each ppm: (We need multiple spectra in order to get a reliable median!)
  m <- matrixStats::colMedians(spectra2)
  # Divide at each ppm by its median:
  f <- spectra2/m[col(spectra2)]
  f <- matrixStats::rowMedians(f)
  # Divide each spectra by its f value
  spectra3 <- spectra2 / f
  list(spectra = spectra3,
       norm_factor = f*areas)
}


#' Normalize NMR samples
#'
#' @param samples A [nmr_dataset] object
#' @param method The criteria to be used for normalization
#' @param values If `method == "value"` then values is a list
#'               with the normalization values. The list must be named as the
#'               data fields to normalize. Typically would be something like:
#'               `values = list(data_1r = c(val1, val2, val3))`.
#'               If `method == "region"` then values is the chemical shift
#'               region to integrate.
#' @return The [nmr_dataset] object, with the samples normalized
#' @export
nmr_normalize <- function(samples,
                          method = c("area", "max", "value", "region", "pqn", "none"),
                          values = NULL) {
  # This function does not consider >1D samples. Some things may work by chance,
  # but it needs testing and revision.
  
  method <- tolower(method[1])
  if (!(method %in% c("area", "max", "value", "region", "pqn", "none"))) {
    stop("Unknown method: ", method)
  }
  if (method == "none") {
    return(samples)
  }
  
  if (samples[["processing"]][["normalization"]]) {
    warning("Samples were already normalized. Applying further normalizations")
  }
  
  data_fields <- names(samples)[grepl(pattern = "^data_.*", x = names(samples))]
  if (method == "area") {
    norm_factor <- list()
    for (data_field in data_fields) {
      norm_factor[[data_field]] <- rowSums(samples[[data_field]])
    }
    samples[["processing"]][["normalization_params"]] <- list(method = method,
                                                              norm_factor = norm_factor)
  } else if (method == "max") {
    norm_factor <- list()
    for (data_field in data_fields) {
      norm_factor[[data_field]] <- apply(samples[[data_field]], 1, max)
    }
    samples[["processing"]][["normalization_params"]] <- list(method = method,
                                                              norm_factor = norm_factor)
  } else if (method == "value") {
    norm_factor <- values
    samples[["processing"]][["normalization_params"]] <- list(method = method,
                                                              norm_factor = norm_factor)
  } else if (method == "region") {
    norm_factor <- list()
    for (data_field in data_fields) {
      if (length(samples[["axis"]]) > 1) {
        stop("Region normalization not implemented for dimensionality > 1")
      }
      axis1 <- samples[["axis"]][[1]]
      region_range <- axis1 >= min(values) & axis1 <= max(values)
      norm_factor[[data_field]] <- rowSums(samples[[data_field]][,region_range])
    }
    samples[["processing"]][["normalization_params"]] <- list(method = method,
                                                              norm_factor = norm_factor)
  } else if (method == "pqn") {
    # only for 1D
    if (length(samples[["axis"]]) > 1) {
      stop("PQN normalization not implemented for dimensionality > 1")
    }
    norm_factor <- list()
    for (data_field in data_fields) {
      norm_result <- norm_pqn(samples[[data_field]])
      samples[[data_field]] <- norm_result$spectra
      norm_factor[[data_field]] <- norm_result$norm_factor
    }
    samples[["processing"]][["normalization"]] <- TRUE
    return(samples)
  } else {
    stop("Unimplemented method: ", method)
  }
  for (data_field in data_fields) {
    samples[[data_field]] <- sweep(x = samples[[data_field]],
                                   MARGIN = 1,
                                   STATS = norm_factor[[data_field]],
                                   FUN = "/")
  }
  samples[["processing"]][["normalization"]] <- TRUE
  return(samples)
}