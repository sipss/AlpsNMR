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
#'
#' The integration is very naïve and consists of the sum of the intensities
#' in the given ppm range.
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
nmr_integrate_regions <- function(samples, regions, fix_baseline = TRUE) {
  UseMethod("nmr_integrate_regions")
}

#' @rdname nmr_integrate_regions
#' @export
nmr_integrate_regions.nmr_dataset <- function(samples, regions, fix_baseline = TRUE) {
  areas <- purrr::map_dfc(regions, function(region) {
    to_sum <- samples$axis[[1]] >= min(region) & samples$axis[[1]] < max(region)
    region_to_sum <- samples$data_1r[, to_sum]
    area <- rowSums(region_to_sum)
    if (fix_baseline) {
      basel <- t(apply(region_to_sum, 1, rough_baseline))
      area_basel <- rowSums(basel)
      area <- area - area_basel
    }
    area
  })
  dplyr::bind_cols(nmr_get_metadata(samples, "NMRExperiment"),
                   areas)
}

rough_baseline <- function(x) {
  n <- length(x)
  basel <- signal::interp1(x = c(1, n),
                           y = x[c(1, n)],
                           xi = seq_len(n))
  basel <- ifelse(basel > x, x, basel)
  basel
}

#' @rdname nmr_integrate_regions
#' @export
nmr_integrate_regions.nmr_dataset_1D <- function(samples, regions, fix_baseline = TRUE) {
  if (is.null(names(regions))) {
    names(regions) <- purrr::map_chr(regions, ~sprintf("ppm_%4.4f", mean(.)))
  }
  areas <- purrr::map_dfc(regions, function(region) {
    to_sum <- samples$axis >= min(region) & samples$axis < max(region)
    region_to_sum <- samples$data_1r[, to_sum]
    area <- rowSums(region_to_sum)
    if (fix_baseline) {
      basel <- t(apply(region_to_sum, 1, rough_baseline))
      area_basel <- rowSums(basel)
      area <- area - area_basel
    }
    area
  })
  dplyr::bind_cols(nmr_get_metadata(samples, "NMRExperiment"),
                   areas)
}