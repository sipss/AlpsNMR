
#' Given a list of inflection points, find the ones closer to a peak apex.
#'
#' @param peak_limit_left The indices of all possible left inflection points
#' @param peak_limit_right  The indices of all possible right inflection points
#' @param pos The index of the apex we are looking at
#'
#' @return a named numeric vector with the indices for the left inflection point, the apex and the right inflection point.
#' @noRd
get_peak_bounds <- function(peak_limit_left, peak_limit_right, pos){
    peak_limit_left <- peak_limit_left[pos - peak_limit_left > 0]
    peak_limit_right <- peak_limit_right[peak_limit_right - pos > 0]
    if(length(peak_limit_right) == 0 || length(peak_limit_left) == 0){
        return(c("left" = 0, "apex" = pos, "right" = 0))
    }
    right <- peak_limit_right[1]
    left <- peak_limit_left[length(peak_limit_left)]
    return(c("left" = left, "apex" = pos, "right" = right))
}


#' A lorentzian function estimation
#' 
#' \eqn{f(x, x_0, \gamma, A) = A \frac{1}{\pi \gamma} \frac{\gamma^2}{(x-x_0)^2 + \gamma^2}}
#'
#' @param x Numbers where the function will be evaluated.
#' @param x0 The peak location
#' @param gamma Lorentzian peak width
#' @param A A multiplication factor. When \eqn{x = x_0}, the lorentzian will be \eqn{\frac{A}{\pi \gamma}}
#'
#' @return A numeric vector of the same length than `x`, with the values of the lorentzian
#' @noRd
lorentzian <- function(x, x0, gamma, A) {
    A * (1/(pi*gamma)) * ((gamma^2) / ((x-x0)^2 + gamma^2))
}


get_norm_rmse <- function(y_fitted, y, y_fitted_apex, y_apex){
    y_fitted_adj <- y_fitted - y_fitted_apex + y_apex
    mse <- sum((y - y_fitted_adj)**2)/length(y)
    rmse <- sqrt(mse)
    norm_rmse <- rmse/y_fitted_apex
    return(norm_rmse)
}


fit_lorentzians_to_peak_data <- function(peak_data, nmr_dataset, nrmse_max = Inf, area_min = 0, area_max = Inf, ppm_min = 0.1, ppm_max = 12) {
    peak_data$gamma <- NA_real_
    peak_data$area <- NA_real_
    peak_data$norm_rmse <- NA_real_
    peak_data$accepted <- NA
    
    
    x <- nmr_dataset$axis
    sindex_prev <- -1
    # For each peak:
    for (i in seq_len(nrow(peak_data))) {
        posi <- peak_data$pos[i]
        sindex <- peak_data$sample_idx[i]
        # We want to estimate a lorentzian in a computationally fast way.
        # peak_data is sorted by sample, so we avoid recalculating.
        if (sindex != sindex_prev) {
            y <- nmr_dataset$data_1r[sindex,]
            sgf <- signal::sgolayfilt(y,p=2, n=21, m = 2, ts=x[2]-x[1]) # 2nd derivative
            peak_limit_left <- which((diff(sign(sgf)))== -2)
            peak_limit_right <- which((diff(sign(sgf)))== 2)
        }
        peak_bounds <- get_peak_bounds(peak_limit_left, peak_limit_right, posi)
        # Estimate gamma with the inflection points:
        # The lorentzian second derivative:
        # $$f''(x, x_0, A, \gamma) = -\frac{2A \gamma \left(\gamma^2-3\left(x-x_0\right)^2\right)}{\pi \left(\gamma^2+\left(x-x_0\right)^2\right)^3}$$
        # The inflection points are at:
        # f''(x, x_0, A, \gamma) = 0 \rightarrow
        #    x_{right}=x_0+\frac{\sqrt{3}\gamma}{3}
        #    x_{left}=x_0-\frac{\sqrt{3}\gamma}{3}
        # Therefore gamma is:
        # \gamma = \frac{\sqrt{3}}{2} (x_{right} - x_{left} )
        # FIXME: Possible minor improvement. Use a linear interpolation to
        #        determine x[peak_bounds["right"]] and the corresponding left.
        #        The error in in gamma would be reduced by half.
        gamma <- sqrt(3)/2 * (x[peak_bounds["right"]] - x[peak_bounds["left"]])
        # Estimate the amplitude of the lorentzian with the amplitude of the second derivative. 
        # This has the effect of removing the baseline up to linear effects.
        # At $x = x_0$, the second derivative is:
        # f''(x = x_0, x_0, A, \gamma) = -\frac{2A}{\pi \gamma^3}
        # So we can estimate A as follows:
        estimated_A <- - (pi * gamma^3 * sgf[posi])/2
        # And then get the whole lorentzian:
        y_fitted <- lorentzian(
            x = x[peak_bounds["left"]:peak_bounds["right"]],
            x0 = x[peak_bounds["apex"]],
            gamma = gamma,
            A = estimated_A
        )
        # So we can estimate how good our fit is to the real data:
        y_fitted_apex_idx <- peak_bounds["apex"] - peak_bounds["left"] + 1L
        norm_rmse <- get_norm_rmse(
            y_fitted,
            y[peak_bounds["left"]:peak_bounds["right"]],
            y_fitted_apex = y_fitted[y_fitted_apex_idx],
            y_fitted = y[peak_bounds["apex"]]
        )
        # And save the result in the peak_data table:
        peak_data$gamma[i] <- gamma
        peak_data$area[i] <- estimated_A
        peak_data$norm_rmse[i] <- norm_rmse
        peak_data$accepted[i] <- (
            norm_rmse <= nrmse_max &&
                estimated_A >= area_min &&
                estimated_A <= area_max &&
                x[peak_bounds["apex"]] >= ppm_min &&
                x[peak_bounds["apex"]] <= ppm_max &&
                !is_ppm_region_excluded(
                    region = c(x[peak_bounds["left"]], x[peak_bounds["right"]]),
                    nmr_get_excluded_regions(nmr_dataset)
                )
        )
        sindex_prev <- sindex
    }
    return(peak_data)
}

