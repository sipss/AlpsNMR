#' PQN normalization
#' @noRd
#' @param spectra A matrix with one spectrum on each row
#' @return A list with `spectra` (a matrix with one spectrum on each row) and `norm_factor`
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
  UseMethod("nmr_normalize")
}

#' @rdname nmr_normalize
#' @family nmr_dataset_1D functions
#' @export
nmr_normalize.nmr_dataset_1D <- function(samples,
                                         method = c("area", "max", "value", "region", "pqn", "none"),
                                         values = NULL) {
  method <- tolower(method[1])
  if (!(method %in% c("area", "max", "value", "region", "pqn", "none"))) {
    stop("Unknown method: ", method)
  }
  if (method == "none") {
    return(samples)
  }
  
  if (method == "area") {
    norm_factor <- rowSums(samples[["data_1r"]])
  } else if (method == "max") {
    norm_factor <- apply(samples[["data_1r"]], 1, max)
  } else if (method == "value") {
    norm_factor <- values
  } else if (method == "region") {
    region_range <- samples[["axis"]] >= min(values) & samples[["axis"]] <= max(values)
    norm_factor <- rowSums(samples[["data_1r"]][,region_range])
  } else if (method == "pqn") {
    norm_result <- norm_pqn(samples[["data_1r"]])
    samples[["data_1r"]] <- norm_result$spectra
    norm_factor <- norm_result$norm_factor
    nmr_diagnose(samples) <- list(name = "Normalization",
                                  norm_factor = norm_factor)
    return(samples)
  } else {
    stop("Unimplemented method: ", method)
  }
  samples[["data_1r"]] <- sweep(x = samples[["data_1r"]],
                                MARGIN = 1,
                                STATS = norm_factor,
                                FUN = "/")
  nmr_diagnose(samples) <- list(name = "Normalization",
                                norm_factor = norm_factor)
  return(samples)
}
