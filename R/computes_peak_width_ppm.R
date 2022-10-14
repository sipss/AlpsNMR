#' Peak width estimation for integration
#'
#' Estimates the peak width (ppm width) to perform peak integration using
#' `nmr_integrate_peak_positions`. for this purpose, the full width at half maximum
#' of a peak from alanine doublet is considered.
#'
#' @family peak integration functions
#' @param nmr_dataset An [nmr_dataset_1D].
#' @return Numerical. A peak width (ppm) that may be set in `nmr_integrate_peak_positions`
#' @noRd
computes_peak_width_ppm <- function(nmr_dataset) {
    ppms <- nmr_dataset$axis
    alanine <- which((nmr_dataset$axis > 1.485) &
        (nmr_dataset$axis < 1.515))
    alanine_vector <- ppms[alanine]

    x <- alanine_vector
    y <- nmr_dataset$data_1r[1, alanine]
    d <- data.frame(x, y)
    # plot(d, type="l", xlab = "Alanine", ylab = "a.u.")
    xmax <- d$x[d$y == max(d$y)]

    # x1 <- d$x[d$x < xmax][which.min(abs(d$y[d$x < xmax]-max(d$y)/2))]
    x2 <- d$x[d$x > xmax][which.min(abs(d$y[d$x > xmax] - max(d$y) / 2))]
    x3 <- abs((x2 - xmax) * 2)

    # points(c(x2, x4), c(d$y[d$x==x2], d$y[d$x==x4]), col="red")
    FWHM <- x3
    message("calculated width for integration is ", FWHM, " ppm")
    return(FWHM)
}
