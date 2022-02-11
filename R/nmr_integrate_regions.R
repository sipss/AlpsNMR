#' Integrate regions
#'
#' Integrate given regions and return a data frame with them
#' @param samples A [nmr_dataset] object
#' @param regions A named list. Each element of the list is a region,
#'                                given as a named numeric vector of length two with the range
#'                                to integrate. The name of the region will be the name of the
#'                                column
#'
#' @return An [nmr_dataset_peak_table] object
#'
#' @examples
#' #Creating a dataset
#' dataset <- new_nmr_dataset_1D(ppm_axis = 1:10,
#'                               data_1r = matrix(sample(0:99,replace = TRUE), nrow = 10),
#'                               metadata = list(external = data.frame(NMRExperiment = c("10", 
#'                               "20", "30", "40", "50", "60", "70", "80", "90", "100"))))
#'
#' # Integrating selected regions
#' peak_table_integration = nmr_integrate_regions(
#'                                    samples = dataset,
#'                                    regions = list(ppm = c(2,5)),
#'                                    fix_baseline = TRUE)
#'
#' @export
#' @family peak detection functions
#' @family peak integration functions
nmr_integrate_regions <- function(samples, regions, ...) {
    UseMethod("nmr_integrate_regions")
}

rough_baseline <- function(x, allow_baseline_above_signal = TRUE) {
    n <- length(x)
    if (n == 0) {
        return(numeric(0L))
    }
    basel <- signal::interp1(x = c(1, n),
                             y = x[c(1, n)],
                             xi = seq_len(n))
    if (!allow_baseline_above_signal) {
        basel <- ifelse(basel > x, x, basel)
    }
    basel
}

#' @rdname nmr_integrate_regions
#' @family nmr_dataset_1D functions
#' @param fix_baseline A logical. If `TRUE` it removes the baseline. See details
#'                                below
#' @param excluded_regions_as_zero A logical. It determines the behaviour of the
#'    integration when integrating regions that have been excluded. If `TRUE`,
#'    it will treat those regions as zero. If `FALSE` (the default) it will return
#'    NA values.
#'
#'    If `fix_baseline` is `TRUE`, then the region boundaries are used to estimate
#'    a baseline. The baseline is estimated "connecting the boundaries with a straight
#'    line". Only when the spectrum is above the baseline the area is integrated
#'    (negative contributions due to the baseline estimation are ignored).
#'
#' @param set_negative_areas_to_zero A logical. Ignored if `fix_baseline` is `FALSE`.
#'     When set to `TRUE` negative areas are set to zero.
#'
#' @param ... Keep for compatibility
#' @export
#' @examples
#' #Creating a dataset
#' dataset <- new_nmr_dataset_1D(ppm_axis = 1:10,
#'                               data_1r = matrix(sample(0:99,replace = TRUE), nrow = 10),
#'                               metadata = list(external = data.frame(NMRExperiment = c("10",
#'                                "20", "30", "40", "50", "60", "70", "80", "90", "100"))))
#' 
#' # Integrating selected regions
#' peak_table_integration = nmr_integrate_regions(
#'                                    samples = dataset,
#'                                    regions = list(ppm = c(2,5)),
#'                                    fix_baseline = TRUE)
#'         
nmr_integrate_regions.nmr_dataset_1D <- function(samples,
                                                 regions,
                                                 fix_baseline = TRUE,
                                                 excluded_regions_as_zero = FALSE,
                                                 set_negative_areas_to_zero = FALSE,
                                                 ...) {
    if (is.null(names(regions))) {
        names(regions) <-
            purrr::map_chr(regions, ~ sprintf("ppm_%4.4f", mean(.)))
    }
    ppm_res <- nmr_ppm_resolution(samples)
    areas <- purrr::map_dfc(regions, function(region) {
        to_sum <- samples$axis >= min(region) & samples$axis < max(region)
        # If there are no ppm to sum, or if we are integrating an excluded region,
        # then return NA.
        if (isTRUE(excluded_regions_as_zero)) {
            integrating_allowed <- TRUE
        } else {
            # I don't care about 3*ppm_res or 2*ppm_res, but I leave some margin just in case:
            integrating_allowed <-
                all(diff(samples$axis[to_sum]) < 3 * ppm_res)
        }
        if (all(to_sum == FALSE) || !integrating_allowed) {
            return(NA * numeric(nrow(samples$data_1r)))
        }
        region_to_sum <- samples$data_1r[, to_sum]
        area <- rowSums(region_to_sum)
        if (fix_baseline) {
            basel <- t(
                apply(
                    region_to_sum,
                    1,
                    rough_baseline,
                    allow_baseline_above_signal = !set_negative_areas_to_zero
                )
            )
            area_basel <- rowSums(basel)
            area <- area - area_basel
        }
        area * ppm_res
    })
    new_nmr_dataset_peak_table(peak_table = as.matrix(areas),
                               metadata = samples$metadata)
}

#' Integrate peak positions
#' 
#' The function allows the integration of a given ppm vector with a specific width.
#' 
#' @return Integrate peak positions
#' @name nmr_integrate_peak_positions
#' @param samples A [nmr_dataset] object
#' @inheritParams regions_from_peak_table
#' @inheritDotParams nmr_integrate_regions
#'
#' @inherit nmr_integrate_regions return
#' @export
#' @family peak integration functions
#' @family nmr_dataset_1D functions
nmr_integrate_peak_positions <- function(samples,
                                         peak_pos_ppm,
                                         peak_width_ppm = 0.006,
                                         ...) {
    # Computes the alanine peak_width_ppm
    if (is.null(peak_width_ppm)) {
        peak_width_ppm <- computes_peak_width_ppm(samples)
    }
    
    # dataframe as input
    if (is.data.frame(peak_pos_ppm)) {
        message("peak_pos_ppm input introduced as dataframe")
        peak_pos_ppm = peak_pos_ppm$ppm
    }
    
    regions <- regions_from_peak_table(peak_pos_ppm, peak_width_ppm)
    nmr_integrate_regions(samples, regions, ...)
}

#' Get integrals with metadata from `integrate peak positions`
#' 
#' @param integration_object A [nmr_dataset] object
#' @return Get integrals with metadata from `integrate peak positions`
#' @export
#'
#' @family peak integration functions
#' @family nmr_dataset_1D functions
#' @return integration dataframe
get_integration_with_metadata <- function(integration_object) {
    integration_data <- nmr_data(integration_object)
    meta_data <- nmr_meta_get(integration_object, groups = "external")
    integration_dataframe <- cbind(meta_data, integration_data)
    return(integration_dataframe)
}
