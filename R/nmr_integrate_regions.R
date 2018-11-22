#' Integrate regions
#' 
#' Integrate given regions and return a data frame with them
#' @param samples A [nmr_dataset] object
#' @param regions A named list. Each element of the list is a region,
#'                given as a named numeric vector of length two with the range
#'                to integrate. The name of the region will be the name of the
#'                column
#' @param fix_baseline A logical. If `TRUE` it removes the baseline. See details
#'                below
#' @param excluded_regions_as_zero A logical. It determines the behaviour of the
#'  integration when integrating regions that have been excluded. If `TRUE`, 
#'  it will treat those regions as zero. If `FALSE` (the default) it will return
#'  NA values.
#'
#' If `fix_baseline` is `TRUE`, then the region boundaries are used to estimate
#' a baseline. The baseline is estimated "connecting the boundaries with a straight
#' line". Only when the spectrum is above the baseline the area is integrated
#' (negative contributions due to the baseline estimation are ignored).
#' 
#' @return A data frame with the NMRExperiment column and one additional column
#'         for each given region.
#'         
#' @examples 
#' \dontrun{
#' # We integrate a region with two peaks and a valley. This is how the
#' # final area is integrated.
#' library(ggplot2)
#' library(signal)
#' x <- seq(from = 1, to = 11, by = 0.05)
#' y <- signal::interp1(x = 1:11,
#'                      y = c(10, 10, 10, 12, 14, 20, 9, 14, 11, 11, 12),
#'                      xi = x)
#' 
#' xb <- c(which.min(abs(x - 3)), which.min(abs(x - 9)))
#' basel <- NIHSnmr:::rough_baseline(y[xb[1]:xb[2]])
#' 
#' ggplot(mapping = aes(x = x, y = y)) +
#'   geom_line(data = data.frame(x = x, y = y)) + 
#'   geom_polygon(data = data.frame(x = x[c(xb[1]:xb[2], rev(xb[1]:xb[2]))],
#'                                  y = c(y[xb[1]:xb[2]], rev(basel))),
#'                fill = "blue") +
#'   scale_y_continuous(limits = c(5, 20))
#' }
#' @export
nmr_integrate_regions <- function(samples, regions, fix_baseline = TRUE,
                                  excluded_regions_as_zero = FALSE) {
  UseMethod("nmr_integrate_regions")
}

rough_baseline <- function(x) {
  n <- length(x)
  if (n == 0) {
    return(numeric(0L))
  }
  basel <- signal::interp1(x = c(1, n),
                           y = x[c(1, n)],
                           xi = seq_len(n))
  basel <- ifelse(basel > x, x, basel)
  basel
}

#' @rdname nmr_integrate_regions
#' @export
nmr_integrate_regions.nmr_dataset_1D <- function(samples, regions, fix_baseline = TRUE,
                                                 excluded_regions_as_zero = FALSE) {
  if (is.null(names(regions))) {
    names(regions) <- purrr::map_chr(regions, ~sprintf("ppm_%4.4f", mean(.)))
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
      integrating_allowed <- all(diff(samples$axis[to_sum]) < 3*ppm_res)
    }
    if (all(to_sum == FALSE) || !integrating_allowed) {
      return(NA*numeric(nrow(samples$data_1r)))
    }
    region_to_sum <- samples$data_1r[, to_sum]
    area <- rowSums(region_to_sum)
    if (fix_baseline) {
      basel <- t(apply(region_to_sum, 1, rough_baseline))
      area_basel <- rowSums(basel)
      area <- area - area_basel
    }
    area*ppm_res
  })
  dplyr::bind_cols(nmr_meta_get(samples, "NMRExperiment"),
                   areas)
}

#' Integrate peak positions
#'
#' @inheritParams nmr_integrate_regions
#' @inheritParams regions_from_peak_table
#'
#' @return A data frame with the NMRExperiment column and one additional column
#'         for each peak position.
#' @export
#'
nmr_integrate_peak_positions <- function(samples,
                                         peak_pos_ppm,
                                         peak_width_ppm,
                                         fix_baseline = TRUE,
                                         excluded_regions_as_zero = FALSE) {
  regions <- regions_from_peak_table(peak_pos_ppm, peak_width_ppm)
  nmr_integrate_regions(samples, regions, fix_baseline = fix_baseline,
                        excluded_regions_as_zero = excluded_regions_as_zero)
}
