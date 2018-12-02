#' Integrate regions
#' 
#' Integrate given regions and return a data frame with them
#' @param samples A [nmr_dataset] object
#' @param regions A named list. Each element of the list is a region,
#'                given as a named numeric vector of length two with the range
#'                to integrate. The name of the region will be the name of the
#'                column
#' 
#' @return An [nmr_dataset_peak_table] object
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
#'                below
#' @param excluded_regions_as_zero A logical. It determines the behaviour of the
#'  integration when integrating regions that have been excluded. If `TRUE`, 
#'  it will treat those regions as zero. If `FALSE` (the default) it will return
#'  NA values.
#'
#'  If `fix_baseline` is `TRUE`, then the region boundaries are used to estimate
#'  a baseline. The baseline is estimated "connecting the boundaries with a straight
#'  line". Only when the spectrum is above the baseline the area is integrated
#'  (negative contributions due to the baseline estimation are ignored).
#' 
#' @param set_negative_areas_to_zero A logical. Ignored if `fix_baseline` is `FALSE`.
#'   When set to `TRUE` negative areas are set to zero.
#'
#' @param ... Keep for compatibility
#' @export
nmr_integrate_regions.nmr_dataset_1D <- function(samples, regions,
                                                 fix_baseline = TRUE,
                                                 excluded_regions_as_zero = FALSE,
                                                 set_negative_areas_to_zero = FALSE,
                                                 ...) {
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
      basel <- t(apply(region_to_sum, 1, 
                       rough_baseline,
                       allow_baseline_above_signal = !set_negative_areas_to_zero))
      area_basel <- rowSums(basel)
      area <- area - area_basel
    }
    area*ppm_res
  })
  new_nmr_dataset_peak_table(peak_table = as.matrix(areas),
                             metadata = samples$metadata)
}

#' Integrate peak positions
#'
#' @param samples A [nmr_dataset] object
#' @inheritParams regions_from_peak_table
#' @inheritDotParams nmr_integrate_regions
#'
#' @inherit nmr_integrate_regions return
#' @export
#'
#' @family peak integration functions
#' @family nmr_dataset_1D functions
nmr_integrate_peak_positions <- function(samples,
                                         peak_pos_ppm,
                                         peak_width_ppm,
                                         ...) {
  regions <- regions_from_peak_table(peak_pos_ppm, peak_width_ppm)
  nmr_integrate_regions(samples, regions, ...)
}
