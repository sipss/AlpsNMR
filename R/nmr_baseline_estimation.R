#' Baseline Removal NMR
#'
#' Removes the baseline on an [nmr_dataset_1D] object, using [baseline::baseline.als].
#'
#' **Deprecated**: `nmr_baseline_removal()` will be removed in early 2027.
#' Use [nmr_baseline_estimation()] instead: it estimates the baseline without
#' overwriting `data_1r`, and downstream functions (e.g.
#' [nmr_baseline_threshold()], [nmr_detect_peaks()]) pick it up automatically
#' when present.
#'
#' @family baseline removal functions
#' @seealso [baseline::baseline.als]
#' @param nmr_dataset An [nmr_dataset_1D].
#' @inheritParams baseline::baseline.als
#' @return The same [nmr_dataset_1D] object after baseline removal.
#' @export
#'
#' @examples
#' dataset_1D <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' dataset_no_base_line <- nmr_baseline_removal(dataset_1D, lambda = 6, p = 0.01)
#'
nmr_baseline_removal <- function(nmr_dataset,
    lambda = 6,
    p = 0.05,
    maxit = 20) {
    cli::cli_warn(
        c(
            "!" = "{.fn nmr_baseline_removal} is deprecated and will be removed in early 2027.",
            "i" = "Use {.fn nmr_baseline_estimation} instead."
        ),
        .frequency = "regularly",
        .frequency_id = "nmr_baseline_removal_deprecated",
    )
    results <- baseline::baseline(
        nmr_dataset$data_1r,
        method = "als",
        lambda = lambda,
        p = p,
        maxit = maxit
    )
    nmr_dataset$data_1r <- baseline::getCorrected(results)
    nmr_dataset
}


#' Estimate the baseline on an [nmr_dataset_1D] object, using [psalsa()]
#'
#' Estimates the baseline of every sample in `nmr_dataset` with the PSALSA
#' algorithm (see [psalsa()]) and stores it in the `data_1r_baseline`
#' element, leaving `data_1r` itself untouched. Several other functions
#' ([nmr_baseline_threshold()], [nmr_detect_peaks()],
#' [nmr_integrate_regions()], [nmr_normalize()]) pick up `data_1r_baseline`
#' automatically when it is present.
#'
#' `lambda`, `p` and `k` each default to `"auto"`. Whenever any of them is
#' `"auto"`, [tune_psalsa()] is run once on every sample in `nmr_dataset`
#' (pooled together, as it would be for a single call to [tune_psalsa()] with
#' a list of spectra) to pick values for every `"auto"` parameter; a
#' parameter given as an explicit number instead bypasses tuning for that
#' parameter and is passed to [psalsa()] as-is. The same `lambda`/`p`/`k` are
#' then used for every sample. `maxit` also defaults to `"auto"`, meaning
#' [psalsa()]'s own default is used, since `maxit` is not tuned by
#' [tune_psalsa()].
#'
#' @family baseline removal functions
#' @seealso [psalsa()], [tune_psalsa()]
#' @param nmr_dataset An [nmr_dataset_1D].
#' @param lambda Smoothing parameter, or `"auto"` to pick it with
#'   [tune_psalsa()]. See [psalsa()].
#' @param p Asymmetry parameter, or `"auto"` to pick it with [tune_psalsa()].
#'   See [psalsa()].
#' @param k Peak height parameter, or `"auto"` to pick it with
#'   [tune_psalsa()]. See [psalsa()].
#' @param maxit Maximum number of iterations, or `"auto"` to use [psalsa()]'s
#'   own default.
#' @return The same [nmr_dataset_1D] object with the `data_1r_baseline` element.
#' @export
#'
#' @examples
#' dataset_1D <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' dataset_1D <- nmr_baseline_estimation(dataset_1D)
#'
nmr_baseline_estimation <- function(nmr_dataset,
    lambda = "auto",
    p = "auto",
    k = "auto",
    maxit = "auto") {
    spectra <- nmr_dataset$data_1r

    if (identical(lambda, "auto") || identical(p, "auto") || identical(k, "auto")) {
        y_list <- lapply(seq_len(nrow(spectra)), function(i) spectra[i, ])
        tuned <- tune_psalsa(y_list)
        if (identical(lambda, "auto")) {
            lambda <- tuned$lambda
        }
        if (identical(p, "auto")) {
            p <- tuned$p
        }
        if (identical(k, "auto")) {
            k <- tuned$k
        }
    }

    psalsa_args <- list(spectra = spectra, lambda = lambda, p = p, k = k)
    if (!identical(maxit, "auto")) {
        psalsa_args$maxit <- maxit
    }
    result <- do.call(psalsa, psalsa_args)

    data_1r_baseline <- result$baseline
    attr(data_1r_baseline, "psalsa_params") <- list(
        lambda = lambda,
        p = p,
        k = k
    )
    nmr_dataset$data_1r_baseline <- data_1r_baseline
    nmr_dataset
}
