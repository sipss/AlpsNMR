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


#' Normalize nmr_dataset_1D samples
#' 
#' The `nmr_normalize` function is used to normalize all the samples according
#' to a given criteria.
#' 
#' The `nmr_normalize_extra_info` function is used to extract additional information
#' after the normalization. Typically, we want to know what was the actual normalization
#' factor applied to each sample. The extra information includes a plot, representing
#' the dispersion of the normalization factor for each sample.
#'
#' @param samples A [nmr_dataset_1D] object
#'                             
#' @param method The criteria to be used for normalization
#'     - area: Normalize to the total area
#'     - max: Normalize to the maximum intensity
#'     - value: Normalize each sample to a user defined value
#'     - region: Integrate a region and normalize each sample to that region
#'     - pqn: Use Probabalistic Quotient Normalization for normalization
#'     - none: Do not normalize at all
#' @param ... Method dependent arguments:
#'     - `method == "value"`:
#'             - `value`: A numeric vector with the normalization values to use
#'     - `method == "region"`:
#'             - `ppm_range`: A chemical shift region to integrate
#'             - `...`: Other arguments passed on to [nmr_integrate_regions]
#' @return The [nmr_dataset_1D] object, with the samples normalized.
#' Further information for diagnostic of the normalization process is also saved
#' and can be extracted by calling `nmr_normalize_extra_info()` afterwards.
#' @family nmr_dataset_1D functions
#' @export
#' @examples
#' nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' nmr_dataset <- nmr_normalize(nmr_dataset, method = "area")
#' norm_dataset <- nmr_normalize(nmr_dataset)
#' norm_dataset$plot
nmr_normalize <- function(samples, 
                          method = c("area", "max", "value", "region", "pqn", "none"),
                          ...) {
    samples <- validate_nmr_dataset_1D(samples)

    method <- tolower(method[1])
    if (!(method %in% c("area", "max", "value", "region", "pqn", "none"))) {
        stop("Unknown method: ", method)
    }
    if (method == "none") {
        return(samples)
    }
    dots <- list(...)
    
    if (method == "area") {
        norm_factor <- rowSums(samples[["data_1r"]])
    } else if (method == "max") {
        norm_factor <- apply(samples[["data_1r"]], 1, max)
    } else if (method == "value") {
        norm_factor <- dots[["values"]]
    } else if (method == "region") {
        ppm_range <- dots[["ppm_range"]]
        norm_factor <- nmr_integrate_regions(samples, regions = list(ic = ppm_range), ...)
        norm_factor <- norm_factor$peak_table[, "ic"]
    } else if (method == "pqn") {
        norm_result <- norm_pqn(samples[["data_1r"]])
        samples[["data_1r"]] <- norm_result$spectra
        norm_factor <- norm_result$norm_factor
        samples <- nmr_normalize_add_extra_info(samples, norm_factor)
        return(samples)
    } else {
        stop("Unimplemented method: ", method)
    }
    samples[["data_1r"]] <- sweep(x = samples[["data_1r"]],
                                  MARGIN = 1,
                                  STATS = norm_factor,
                                  FUN = "/")
    samples <- nmr_normalize_add_extra_info(samples, norm_factor)
    return(samples)
}

#' @rdname nmr_normalize
#' @export
#' @examples
#' nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' nmr_dataset <- nmr_normalize(nmr_dataset, method = "area")
#' norm_extra_info <- nmr_normalize_extra_info(nmr_dataset)
#' norm_extra_info$plot
nmr_normalize_extra_info <- function(samples) {
    attr(samples, "normalize_extra_info")
}

nmr_normalize_add_extra_info <- function(samples, norm_factor) {
    norm_factor_df <- cbind(nmr_meta_get(samples, columns = "NMRExperiment"),
                                                    norm_factor = norm_factor)
    
    norm_factor_df$norm_factor_norm <- norm_factor_df$norm_factor/stats::median(norm_factor_df$norm_factor)
    
    # This plot should have the y axis in log2 scale
    gplt <- ggplot2::ggplot(norm_factor_df) +
        ggplot2::geom_col(ggplot2::aes(x = .data$NMRExperiment, y = .data$norm_factor_norm)) +
        ggplot2::scale_x_discrete(name = "NMRExperiment", limits = rev(sort(norm_factor_df$NMRExperiment))) +
        ggplot2::geom_hline(yintercept = 1, color = "red") +
        ggplot2::scale_y_continuous(name = "Normalization factors (normalized to the median)") +
        ggplot2::coord_flip()
    attr(samples, "normalize_extra_info") <- list(
        norm_factor = norm_factor_df,
        plot = gplt
    )
    samples
}
