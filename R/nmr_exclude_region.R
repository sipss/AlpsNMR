#' Exclude region from samples
#'
#' Excludes a given region (for instance to remove the water peak)
#' 
#' @param samples An object
#' @param exclude A list with regions to be removed Typically:
#'                                `exclude = list(water = c(4.7, 5.0))`
#' @return The same object, with the regions excluded
#' @export
#' @examples 
#' nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' exclude_regions <- list(water = c(5.1, 4.5))
#' nmr_dataset <- nmr_exclude_region(nmr_dataset, exclude = exclude_regions)
#' 
nmr_exclude_region <- function(samples, exclude = list(water = c(4.7, 5.0))) {
    UseMethod("nmr_exclude_region")
}

#' @rdname nmr_exclude_region
#' @family nmr_dataset_1D functions
#' @export
#' @examples 
#' nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' exclude_regions <- list(water = c(5.1, 4.5))
#' nmr_dataset <- nmr_exclude_region(nmr_dataset, exclude = exclude_regions)
#' 
nmr_exclude_region.nmr_dataset_1D <- function(samples, exclude = list(water = c(4.7, 5.0))) {
    if (is.null(exclude) || length(exclude) == 0) {
        return(samples)
    }
    axis_include <- is_ppm_included(samples[["axis"]], exclude)
    samples[["axis"]] <- samples[["axis"]][axis_include]
    samples[["data_1r"]] <- samples[["data_1r"]][, axis_include, drop=FALSE]
    samples[["excluded_regions"]] <- c(nmr_get_excluded_regions(samples), exclude)
    return(samples)
}

is_ppm_included <- function(ppm, exclude) {
    ppms_included <- rep(TRUE, length(ppm))
    for (i_region in seq_along(exclude)) {
        region <- exclude[[i_region]]
        excl_dim1 <- ppm >= min(region) & ppm <= max(region)
        ppms_included[excl_dim1] <- FALSE
    }
    names(ppms_included) <- ppm
    ppms_included
}

are_ppm_regions_excluded <- function(regions_min, regions_max, exclude) {
    stopifnot(length(regions_min) == length(regions_max))
    stopifnot(regions_min <= regions_max)
    out <- rep(FALSE, length(regions_min))
    for (excl in exclude) {
        excl_min <- min(excl)
        excl_max <- max(excl)
        out[regions_min <= excl_max & regions_max >= excl_min] <- TRUE
    }
    out
}

nmr_get_excluded_regions <- function(samples) {
    samples <- validate_nmr_dataset_1D(samples)
    samples[["excluded_regions"]]
}
