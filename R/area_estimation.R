
#' Given a list of inflection points, find the ones closer to a peak apex.
#'
#' @param peak_limit_left The indices of all possible left inflection points
#' @param peak_limit_right  The indices of all possible right inflection points
#' @param pos The index of the apex we are looking at
#'
#' @return a named numeric vector with the indices for the left inflection point, the apex and the right inflection point.
#' @noRd
get_peak_bounds <- function(peak_limit_left, peak_limit_right, pos, x, sgf){
    peak_limit_left <- peak_limit_left[pos - peak_limit_left > 0]
    peak_limit_right <- peak_limit_right[peak_limit_right - pos > 0]
    if(length(peak_limit_right) == 0 || length(peak_limit_left) == 0) {
        return(c("left" = 0, "apex" = pos, "right" = 0, "xleft" = 0, "xright" = 0))
    }
    right <- peak_limit_right[1]
    left <- peak_limit_left[length(peak_limit_left)]
    # y = y0 + m*(x-x0)
    # y = 0 ==> x = x0 -1/m * y0
    mfun <- function(x0, y0, x1, y1) {(y1-y0)/(x1-x0)}
    
    if (left == length(sgf)) {
        xleft <- x[left]
    } else {
        mleft <- mfun(x[left], sgf[left], x[left+1L], sgf[left+1L])
        xleft <- x[left] - 1/mleft * sgf[left]
    }
    
    if (right == length(sgf)) {
        xright <- x[right]
    } else {
        mright <- mfun(x[right], sgf[right], x[right+1L], sgf[right+1L])
        xright <- x[right] - 1.0/mright * sgf[right]
    }
    return(c("left" = left, "apex" = pos, "right" = right, "xleft" = xleft, "xright" = xright))
}


#' A lorentzian function estimation
#' 
#' \eqn{f(x, x_0, \gamma, A) = A \frac{1}{\pi \gamma} \frac{\gamma^2}{(x-x_0)^2 + \gamma^2}}
#'
#' @param x Numbers where the function will be evaluated.
#' @param x0 The peak location
#' @param gamma Lorentzian peak width (in the same units as x and x0)
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


refine_lorentzian_fit_with_nls <- function(data_to_fit, start, method) {
    estimated_A <- start[["A"]]
    gamma <- start[["gamma"]]
    peak_pos_ppm <- start[["x0"]]
    error_msgs <- c()
    
    if (method == "peak") {
        formula <- y ~ A * (1/(pi*gamma)) * ((gamma^2) / ((x-x0)^2 + gamma^2))
    } else if (method == "2nd_derivative") {
        formula <- y ~ -2*A*gamma/pi * (gamma^2 - 3*(x-x0)^2) / ( gamma^2 + (x-x0)^2)^3
    } else {
        rlang::abort("Unknown method")
    }
    
    tryCatch({
        mm.model.nls <- stats::nls(
            formula = formula,
            data = data_to_fit,
            start = start
        )
        gamma <- mm.model.nls$m$getAllPars()["gamma"]
        estimated_A <- mm.model.nls$m$getAllPars()["A"]
        if (gamma < 0 && estimated_A < 0) {
            gamma <- -gamma
            estimated_A <- -estimated_A
        }
        if (abs(start[["x0"]] - mm.model.nls$m$getAllPars()["x0"]) > gamma) {
            error_msgs <- c(error_msgs, "peak position moved more than gamma!")
        }
        peak_pos_ppm <- mm.model.nls$m$getAllPars()["x0"]
        
    }, error = function(e) {
        msg <- conditionMessage(e)
        error_msgs <- c(error_msgs, msg)
    })
    list(
        estimated_A = estimated_A,
        gamma = gamma,
        peak_pos_ppm = peak_pos_ppm,
        error_msgs = error_msgs
    )
}



#' Fit lorentzians to each peak to estimate areas
#' 
#' The different methods are available for benchmarking while developing, we should pick one.
#' 
#' - gamma is estimated using the inflection points of the signal and fitting them to the lorentzian inflection points
#' - $A$ is estimated using the `amplitude_method` below
#' - The peak position ($x_0$) is given in `peak_data`
#' 
#' Those estimations may be refined with non-linear least squares using `refine_peak_model`. If the nls does not converge,
#' the initial estimations are kept. Convergence -and other nls errors- are saved for further reference and diagnostic.
#' Use `attr(peak_data_fitted, "errors")` to retreive the error messages, where `peak_data_fitted` is assumed to be the
#' output of this function. The refining improves gamma, $A$ and $x_0$.
#' 
#' The baseline estimation (when calculated, see the arguments) is set to Asymmetric Least Squares with 
#' lambda = 6, p=0.05, maxit=20 and it is probably not optimal... yet.
#'
#'
#'
#' @param peak_data The peak data
#' @param nmr_dataset The nmr_dataset object with the data. This function for now assumes nmr_dataset is NOT be baseline corrected
#' @param amplitude_method The method to estimate the amplitude. It may be:
#'  - `"intensity"`. The amplitude of the peak is proportional to the raw intensity at the apex. This is a bad estimation if
#'    the intensity includes a baseline, because the amplitude of the peak will be overestimated
#'  - `"2nd_derivative"`: The amplitude of the peak is proportional to the second derivative of the raw intensity signal at the apex.
#'    This method aims to correct the "intensity" method, since it is expected that the baseline will be mostly removed
#'    when considering the 2nd derivative of the spectrum. The 2nd derivative is calculated with a 2nd order Savitzky-Golay filter of 21 points.
#'  - `"intensity_without_baseline"`: A baseline is estimated on the whole spectra and subtracted from it. Then the peak amplitude
#'    is proportional to the corrected intensity at the apex (as in the "intensity" method).
#' @param refine_peak_model Whether a non linear least squares fitting should be used to refine the estimated parameters. It can be:
#'  - `"none"`: Do not refine using nls.
#'  - `"peak"`: Use a lorentzian peak model and the baseline corrected spectra.
#'  - `"2nd_derivative"`:
#'
#'
#' 
#' @export
#' @return The given data frame `peak_data`, with added columns:
#'  - inflection points, 
#'  - gamma
#'  - area 
#'  - a norm_rmse fitting error
#'  
#'  As well as some attributes
#'  - "errors": A data frame with any error in the peak fitting
#'  - "fit_baseline": Whether the method used has any consideration for the baseline of the signal (maybe not very useful attribute)
#'  - "method_description": A textual description of what we did, to include it in plots
peaklist_fit_lorentzians <- function(
        peak_data, 
        nmr_dataset, 
        amplitude_method = c("intensity", "2nd_derivative", "intensity_without_baseline"),
        refine_peak_model = c("none", "peak", "2nd_derivative")) {
    peak_data$ppm_infl_min <- NA_real_
    peak_data$ppm_infl_max <- NA_real_
    peak_data$gamma_ppb <- NA_real_
    peak_data$area <- NA_real_
    peak_data$norm_rmse <- NA_real_
    amplitude_method <- match.arg(amplitude_method)
    refine_peak_model <- match.arg(refine_peak_model)
    method_description <- c()
    
    x <- nmr_dataset$axis
    sindex_prev <- -1
    all_errors <- list(peak_id = character(0L), error_msg = character(0L))
    pb <- progress_bar_new(name = "Fitting lorentzians", total = nrow(peak_data))

    # For each peak:
    has_warned_baseline <- FALSE
    for (i in seq_len(nrow(peak_data))) {
        progress_bar_update(pb)
        posi <- peak_data$pos[i]
        sindex <- peak_data$sample_idx[i]
        peak_pos_ppm <- peak_data$ppm[i]
        # We want to estimate a lorentzian in a computationally fast way.
        # peak_data is sorted by sample, so we avoid recalculating.
        if (sindex != sindex_prev) {
            y <- nmr_dataset$data_1r[sindex,]
            if (amplitude_method == "intensity_without_baseline" || refine_peak_model == "peak") {
                if ("data_1r_baseline" %in% names(unclass(nmr_dataset))) {
                    y_basel <- as.numeric(nmr_dataset$data_1r_baseline[sindex,])
                } else {
                    if (!has_warned_baseline) {
                        rlang::warn(
                            c(
                                "Estimating the baseline using ALS with lambda 9 and p = 0.05...",
                                "i" = "Use nmr_baseline_estimation before calling this function to customize"
                            )
                        )
                        has_warned_baseline <- TRUE
                    }
                    y_basel <- baseline::baseline.als(
                        spectra = nmr_dataset$data_1r[sindex,, drop = FALSE],
                        lambda = 9,
                        p = 0.05,
                        maxit = 20
                    )$baseline %>% as.numeric()
                }
                y_nobasel <- y - y_basel
            }
            sgf <- signal::sgolayfilt(y, p = 2, n = 21, m = 2, ts = x[2] - x[1]) # 2nd derivative
            peak_limit_left <- which((diff(sign(sgf))) == -2)
            peak_limit_right <- which((diff(sign(sgf))) == 2)
        }
        peak_bounds <- get_peak_bounds(peak_limit_left, peak_limit_right, posi, x, sgf)
        # Estimate gamma with the inflection points:
        # The lorentzian second derivative:
        # $$f''(x, x_0, A, \gamma) = -\frac{2A \gamma \left(\gamma^2-3\left(x-x_0\right)^2\right)}{\pi \left(\gamma^2+\left(x-x_0\right)^2\right)^3}$$
        # The inflection points are at:
        # f''(x, x_0, A, \gamma) = 0 \rightarrow
        #    x_{right}=x_0+\frac{\sqrt{3}\gamma}{3}
        #    x_{left}=x_0-\frac{\sqrt{3}\gamma}{3}
        # Therefore gamma is:
        # \gamma = \frac{\sqrt{3}}{2} (x_{right} - x_{left} )
        gamma <- sqrt(3)/2 * (peak_bounds["xright"] - peak_bounds["xleft"])
        if (amplitude_method == "intensity") {
            # Estimate the amplitude of the lorentzian with the amplitude of the peak 
            estimated_A <- pi*gamma*peak_data$intensity[i]
        }  else if (amplitude_method == "2nd_derivative") {
            # Estimate the amplitude of the lorentzian with the amplitude of the second derivative. 
            # This has the effect of removing the baseline up to linear effects.
            # At $x = x_0$, the second derivative is:
            # f''(x = x_0, x_0, A, \gamma) = -\frac{2A}{\pi \gamma^3}
            # So we can estimate A as follows:
            estimated_A <- -(pi * gamma^3 * sgf[posi])/2
        } else if (amplitude_method == "intensity_without_baseline") {
            # Estimate A based on the amplitude of the signal without baseline
            estimated_A <- y_nobasel[posi]*pi*gamma
        } else {
            rlang::abort(sprintf("amplitude_method '%s'unknown", amplitude_method))
        }
        if (identical(refine_peak_model, "peak")) {
            # Further fitting with nls:
            new_params <- refine_lorentzian_fit_with_nls(
                data_to_fit = data.frame(
                    x = x[peak_bounds["left"]:peak_bounds["right"]],
                    y = y_nobasel[peak_bounds["left"]:peak_bounds["right"]]
                ),
                start = list(
                    A = estimated_A,
                    x0 = x[peak_bounds["apex"]],
                    gamma = gamma
                ),
                method = refine_peak_model
            )
            gamma <- new_params[["gamma"]]
            estimated_A <- new_params[["estimated_A"]]
            peak_pos_ppm <- new_params[["peak_pos_ppm"]]
            if (!is.null(new_params[["error_msgs"]])) {
                all_errors$peak_id <- c(all_errors$peak_id, peak_data$peak_id[i])
                all_errors$error_msg <- paste(new_params[["error_msgs"]], collapse = "\n")
            }
        } else if (identical(refine_peak_model, "2nd_derivative")) {
            # Further fitting with nls:
            new_params <- refine_lorentzian_fit_with_nls(
                data_to_fit = data.frame(
                    x = x[peak_bounds["left"]:peak_bounds["right"]],
                    y = sgf[peak_bounds["left"]:peak_bounds["right"]]
                ),
                start = list(
                    A = estimated_A,
                    x0 = x[peak_bounds["apex"]],
                    gamma = gamma
                ),
                method = refine_peak_model
            )
            gamma <- new_params[["gamma"]]
            estimated_A <- new_params[["estimated_A"]]
            peak_pos_ppm <- new_params[["peak_pos_ppm"]]
            if (!is.null(new_params[["error_msgs"]])) {
                all_errors$peak_id <- c(all_errors$peak_id, peak_data$peak_id[i])
                all_errors$error_msg <- paste(new_params[["error_msgs"]], collapse = "\n")
            }
        } else if (!identical(refine_peak_model, "none")) {
            rlang::abort("Unknown refine_peak_model")
        }
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
        peak_data$ppm_infl_min[i] <- peak_bounds["xleft"]
        peak_data$ppm_infl_max[i] <- peak_bounds["xright"]
        peak_data$gamma_ppb[i] <- 1000*gamma
        peak_data$area[i] <- estimated_A
        peak_data$norm_rmse[i] <- norm_rmse
        sindex_prev <- sindex
    }
    
    # Method description (useful for plotting)
    method_description <- c(method_description, "Estimate \u03b3 with inflection points")
    if (amplitude_method == "intensity") {
        fits_on_baseline <- FALSE
        method_description <- c(
            method_description,
            "Estimate A with raw peak intensity",
            "May have baseline issues"
        )
    }  else if (amplitude_method == "2nd_derivative") {
        fits_on_baseline <- TRUE
        method_description <- c(
            method_description,
            "Estimate A with second derivative",
            "2nd derivative should remove baseline"
        )
    } else if (amplitude_method == "intensity_without_baseline") {
        fits_on_baseline <- TRUE
        method_description <- c(
            method_description,
            "Remove baseline using ALS",
            "Estimate A with the amplitude of the baseline-corrected signal"
        )
        
    }
    
    if (identical(refine_peak_model, "peak")) {
        method_description <- c(
            method_description,
            "Refine \u03b3, A and x0 (nls fits a lorentzian)"
        )
    } else if (identical(refine_peak_model, "2nd_derivative")) {
        method_description <- c(
            method_description,
            "Refine \u03b3, A and x0 (nls fits a lorentzian on 2nd deriv)"
        )
    } else if (identical(refine_peak_model, "none")) {
        method_description <- c(
            method_description,
            "nls refining disabled"
        )
    }
    
    
    
    attr(peak_data, "errors") <- as.data.frame(all_errors)
    attr(peak_data, "fits_on_baseline") <- fits_on_baseline
    attr(peak_data, "method_description") <- paste(
        method_description,
        collapse = "\n"
    )
    peak_data
}




#' Peak list: Create an `accepted` column based on some criteria
#'
#' @param peak_data The peak list (a data frame)
#' @param nmr_dataset The nmr_dataset where the peak_data was computed from
#' @param nrmse_max The normalized root mean squared error of the lorentzian peak fitting must be less than or equal to this value
#' @param area_min Peak areas must be larger or equal to this value
#' @param area_max Peak areas must be smaller or equal to this value
#' @param ppm_min The peak apex must be above this value
#' @param ppm_max The peak apex must be below this value
#' @param keep_rejected If `FALSE`, removes those peaks that do not satisfy the criteria and remove the accepted column (since all would be accepted)
#' @param verbose Print informational message
#'
#' @return The `peak_data`, with a new `accepted` column (or maybe some filtered rows)
#' @export
#'
#' @examples
#' # Fake data:
#' nmr_dataset <- new_nmr_dataset_1D(
#'   1:10,
#'   matrix(c(1:5, 4:2, 3, 0), nrow=1),
#'   list(external = data.frame(NMRExperiment = "10"))
#' )
#' peak_data <- data.frame(
#'   peak_id = c("Peak1", "Peak2"),
#'   NMRExperiment = c("10", "10"),
#'   sample_idx = c(1L,1L),
#'   ppm = c(5, 9),
#'   pos = c(5, 9),
#'   intensity = c(5, 3),
#'   ppm_infl_min = c(3, 8),
#'   ppm_infl_max = c(7, 10),
#'   gamma_ppb = c(1, 1),
#'   area = c(25, 3),
#'   norm_rmse = c(0.01, 0.8)
#' )
#' # Create the accepted column:
#' peak_data <- peaklist_accept_peaks(peak_data, nmr_dataset, area_min = 10, keep_rejected = FALSE)
#' stopifnot(identical(peak_data$peak_id, "Peak1"))
peaklist_accept_peaks <- function(peak_data, nmr_dataset, nrmse_max = Inf, area_min = 0, area_max = Inf, ppm_min = -Inf, ppm_max = Inf, keep_rejected = TRUE, verbose = FALSE) {
    peak_data$accepted <- (
        peak_data$norm_rmse <= nrmse_max &
            peak_data$area >= area_min &
            peak_data$area <= area_max &
            peak_data$ppm >= ppm_min &
            peak_data$ppm <= ppm_max &
            !are_ppm_regions_excluded(peak_data$ppm_infl_min, peak_data$ppm_infl_max, nmr_get_excluded_regions(nmr_dataset))
    )
    report <- c(
        "Acceptance report",
        "i" = glue::glue("{sum(peak_data$accepted)}/{nrow(peak_data)} peaks accepted. ({signif(100*sum(peak_data$accepted)/nrow(peak_data), 3)}%)")
    )
    if (!keep_rejected) {
        report <- c(report, "i" = glue::glue("Removing {sum(!peak_data$accepted)} peaks"))
        peak_data <- peak_data[peak_data$accepted, , drop = FALSE]
        peak_data$accepted <- NULL
    }
    if (verbose) {
        rlang::inform(message = report)
    }
    peak_data
}

