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
#' `"auto"`, `tune_psalsa()` is run once on every sample in `nmr_dataset`
#' (pooled together, as it would be for a single call to `tune_psalsa()` with
#' a list of spectra) to pick values for every `"auto"` parameter; a
#' parameter given as an explicit number (or, when `num_regions` is set, an
#' explicit vector -- see below) instead bypasses tuning for that parameter
#' and is passed to [psalsa()] as-is. The same `lambda`/`p`/`k` are then used
#' for every sample. `maxit` also defaults to `"auto"`, meaning [psalsa()]'s
#' own default is used, since `maxit` is not tuned by `tune_psalsa()`.
#'
#' Whatever `lambda`/`p`/`k` end up being used (auto-tuned or given
#' explicitly), together with the `num_regions` this call was given, are
#' recorded as a `psalsa_params` attribute on the returned `data_1r_baseline`
#' (`attr(dataset$data_1r_baseline, "psalsa_params")`), so `lambda`/`p`/`k`
#' can be reused directly on another, similarly-sized dataset without tuning
#' again: `nmr_baseline_estimation(other_dataset, lambda =
#' psalsa_params$lambda, p = psalsa_params$p, k = psalsa_params$k)`.
#'
#' @family baseline removal functions
#' @seealso [psalsa()], `tune_psalsa()`, `tune_psalsa_spatial()`
#' @param nmr_dataset An [nmr_dataset_1D].
#' @param lambda Smoothing parameter, or `"auto"` to pick it with
#'   `tune_psalsa()`/`tune_psalsa_spatial()`. See [psalsa()].
#' @param p Asymmetry parameter, or `"auto"` to pick it with
#'   `tune_psalsa()`/`tune_psalsa_spatial()`. See [psalsa()].
#' @param k Peak height parameter, or `"auto"` to pick it with
#'   `tune_psalsa()`/`tune_psalsa_spatial()`. See [psalsa()].
#' @param maxit Maximum number of iterations, or `"auto"` to use [psalsa()]'s
#'   own default.
#' @param num_regions If `NULL` (the default), an `"auto"` parameter is tuned
#'   with `tune_psalsa()`, a single signal-wide value shared by every point.
#'   If set, an `"auto"` parameter is instead tuned with
#'   `tune_psalsa_spatial(num_regions = num_regions)`, giving a smooth,
#'   position-varying profile (one value per point) that can differ across
#'   the spectrum -- useful when peak density/baseline behaviour varies
#'   noticeably from region to region. Ignored if `lambda`, `p` and `k` are
#'   all given explicitly (nothing left to tune).
#' @param damping `"auto"` (the default) uses [psalsa()]'s own default
#'   (`damping = 1`, undamped). Set below `1` to under-relax the asymmetric
#'   reweighting update -- see [psalsa()]'s own `damping` parameter for why
#'   this is sometimes necessary: the reweighting depends on a hard
#'   above/below-baseline threshold, and near a flat/quiet region (many
#'   points close to that threshold) the undamped update can oscillate
#'   rather than converge, right up to `maxit`. A `damping` below `1`
#'   typically needs a correspondingly higher `maxit` to still converge.
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
    maxit = "auto",
    num_regions = NULL,
    damping = "auto") {
    spectra <- nmr_dataset$data_1r

    if (identical(lambda, "auto") || identical(p, "auto") || identical(k, "auto")) {
        y_list <- lapply(seq_len(nrow(spectra)), function(i) spectra[i, ])
        tuned <- if (is.null(num_regions)) {
            tune_psalsa(y_list)
        } else {
            tune_psalsa_spatial(y_list, num_regions = num_regions)
        }
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
    if (!identical(damping, "auto")) {
        psalsa_args$damping <- damping
    }
    result <- do.call(psalsa, psalsa_args)

    data_1r_baseline <- result$baseline
    attr(data_1r_baseline, "psalsa_params") <- list(
        lambda = lambda,
        p = p,
        k = k,
        num_regions = num_regions
    )
    nmr_dataset$data_1r_baseline <- data_1r_baseline
    nmr_dataset
}


#' Plot the estimated baseline against the original signal
#'
#' Plots a few samples' original spectra together with their estimated
#' baseline (see [nmr_baseline_estimation()]), so the baseline estimate can
#' be visually inspected against the signal it was estimated from.
#'
#' @family baseline removal functions
#' @seealso [nmr_baseline_estimation()]
#' @param nmr_dataset An [nmr_dataset_1D] with a `data_1r_baseline` element
#'   (i.e. after calling [nmr_baseline_estimation()]).
#' @param NMRExperiment A character vector with the NMRExperiments to plot.
#'   `NULL` (the default) or `"all"` plots every sample; for a legible plot,
#'   pass a handful of NMRExperiments.
#' @param chemshift_range Either a numeric vector of length 2 (a single ppm
#'   range to plot), or a named list of such vectors to plot several regions,
#'   one facet per region, named after the list names.
#' @param nrow,ncol Number of rows/columns of region facets per page. `NULL`
#'   (the default) picks a snug grid for the number of regions requested:
#'   1x`n` for fewer than 4 regions, 2x2 for 4, 2x3 for 5-6, and a fixed 3x3
#'   (paginated via `page`) for 7 or more.
#' @param page Which page of region facets to plot (1-indexed). Requesting a
#'   page beyond the number available is an error.
#' @return A ggplot2 plot with one facet per `chemshift_range` region
#'   (`ggplot2::facet_wrap()`); every requested `NMRExperiment` is overlaid
#'   within each facet. The original signal is drawn as a solid line and the
#'   estimated baseline as a dashed line, both coloured by `NMRExperiment`.
#' @export
#' @examples
#' dataset_1D <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' dataset_1D <- nmr_baseline_estimation(dataset_1D, lambda = 6, p = 0.05)
#' nmr_baseline_estimation_plot(
#'     dataset_1D,
#'     NMRExperiment = names(dataset_1D)[1:2],
#'     chemshift_range = list(Region1 = c(1.2, 1.4), Region2 = c(3.4, 3.6))
#' )
#'
nmr_baseline_estimation_plot <- function(nmr_dataset, NMRExperiment = NULL, chemshift_range = NULL,
    nrow = NULL, ncol = NULL, page = 1) {
    if (!"data_1r_baseline" %in% names(unclass(nmr_dataset))) {
        cli::cli_abort(
            message = c(
                "nmr_dataset has no estimated baseline",
                "i" = "Run {.fn nmr_baseline_estimation} on it first."
            )
        )
    }
    if (is.null(chemshift_range)) {
        cli::cli_abort(
            message = c(
                "chemshift_range must be given",
                "i" = "Pass a numeric vector of length 2 (a single ppm range), or a named list of such vectors for several regions."
            )
        )
    }
    if (is.numeric(chemshift_range)) {
        chemshift_range <- list(chemshift_range)
    }
    if (is.null(names(chemshift_range)) || !all(nzchar(names(chemshift_range)))) {
        names(chemshift_range) <- paste0("Region ", seq_along(chemshift_range))
    }
    axis_range <- range(nmr_dataset$axis)
    for (region_name in names(chemshift_range)) {
        region_i <- chemshift_range[[region_name]]
        if (length(region_i) != 2) {
            cli::cli_abort("Each chemshift_range must have length 2")
        }
        if (max(region_i) < axis_range[1] || min(region_i) > axis_range[2]) {
            cli::cli_abort(
                message = c(
                    "chemshift_range {.val {region_name}} = [{min(region_i)}, {max(region_i)}] ppm falls outside the dataset's axis range",
                    "i" = "The dataset's axis spans [{axis_range[1]}, {axis_range[2]}] ppm."
                )
            )
        }
    }

    if (is.null(NMRExperiment) || identical(NMRExperiment, "all")) {
        NMRExperiment <- names(nmr_dataset)
    }

    num_regions <- length(chemshift_range)
    if (is.null(nrow) && is.null(ncol)) {
        if (num_regions >= 7) {
            nrow <- 3
            ncol <- 3
        } else if (num_regions %in% c(5, 6)) {
            nrow <- 2
            ncol <- 3
        } else if (num_regions == 4) {
            nrow <- 2
            ncol <- 2
        } else {
            nrow <- 1
            ncol <- max(num_regions, 1)
        }
    }
    regions_per_page <- nrow * ncol
    total_pages <- max(ceiling(num_regions / regions_per_page), 1)
    if (page < 1 || page > total_pages) {
        cli::cli_abort(
            message = c(
                "page ({page}) is out of bounds",
                "i" = "There {cli::qty(total_pages)} {?is/are} {total_pages} page{?s} for {num_regions} region{?s} with nrow = {nrow}, ncol = {ncol} ({regions_per_page} per page)."
            )
        )
    }
    page_start <- (page - 1) * regions_per_page + 1
    page_end <- min(page * regions_per_page, num_regions)
    chemshift_range <- chemshift_range[page_start:page_end]

    regions_data <- purrr::imap(chemshift_range, function(range_i, region_name) {
        signal <- tidy(
            nmr_dataset,
            chemshift_range = range_i,
            NMRExperiment = NMRExperiment,
            matrix_name = "data_1r"
        )
        signal$type <- "Signal"
        baseline <- tidy(
            nmr_dataset,
            chemshift_range = range_i,
            NMRExperiment = NMRExperiment,
            matrix_name = "data_1r_baseline"
        )
        baseline$type <- "Baseline"
        region_df <- rbind(signal, baseline)
        region_df$region <- region_name
        region_df
    })
    to_plot <- dplyr::bind_rows(regions_data)
    to_plot$NMRExperiment <- factor(to_plot$NMRExperiment, levels = intersect(NMRExperiment, unique(to_plot$NMRExperiment)))
    to_plot$region <- factor(to_plot$region, levels = names(chemshift_range))
    to_plot$type <- factor(to_plot$type, levels = c("Signal", "Baseline"))

    gplt <- ggplot2::ggplot(
        to_plot,
        ggplot2::aes(
            x = .data$chemshift,
            y = .data$intensity,
            colour = .data$NMRExperiment,
            linetype = .data$type,
            group = interaction(.data$NMRExperiment, .data$type)
        )
    ) +
        ggplot2::geom_line() +
        ggplot2::scale_linetype_manual(values = c(Signal = "solid", Baseline = "dashed"), name = NULL) +
        ggplot2::labs(x = "Chemical Shift (ppm)", y = "Intensity (a.u.)", colour = "NMRExperiment") +
        ggplot2::scale_x_reverse() +
        ggplot2::scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_si(""))) +
        ggplot2::facet_wrap(~ .data$region, nrow = nrow, ncol = ncol, scales = "free")
    gplt
}
