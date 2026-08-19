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
#'   range to plot), or a named list of such vectors to plot several regions
#'   side by side, one facet column per region, named after the list names.
#' @param nrow,ncol Number of NMRExperiments (facet rows) and regions (facet
#'   columns) to show per page. `NULL` (the default) shows every requested
#'   NMRExperiment (`nrow`) and every requested region (`ncol`) on a single
#'   page; set either to paginate that dimension instead of cramming (or
#'   silently subsampling) everything onto one page.
#' @param page Which page to plot (1-indexed). Pages are laid out with
#'   NMRExperiment pages varying fastest, then region pages. Requesting a
#'   page beyond the number available is an error.
#' @return A ggplot2 plot: one facet row per `NMRExperiment` and, when more
#'   than one region is requested, one facet column per `chemshift_range`
#'   region. The original signal is drawn as a solid line and the estimated
#'   baseline as a dashed line, both coloured by `NMRExperiment`.
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

    num_samples <- length(NMRExperiment)
    num_regions <- length(chemshift_range)
    if (is.null(nrow)) {
        nrow <- num_samples
    }
    if (is.null(ncol)) {
        ncol <- num_regions
    }
    num_sample_pages <- max(ceiling(num_samples / nrow), 1)
    num_region_pages <- max(ceiling(num_regions / ncol), 1)
    total_pages <- num_sample_pages * num_region_pages
    if (page < 1 || page > total_pages) {
        cli::cli_abort(
            message = c(
                "page ({page}) is out of bounds",
                "i" = "There {cli::qty(total_pages)} {?is/are} {total_pages} page{?s} for {num_samples} NMRExperiment{?s} and {num_regions} region{?s} with nrow = {nrow}, ncol = {ncol} ({nrow * ncol} facets per page)."
            )
        )
    }
    # NMRExperiment pages vary fastest, then region pages.
    sample_page <- ((page - 1) %% num_sample_pages) + 1
    region_page <- ((page - 1) %/% num_sample_pages) + 1
    sample_start <- (sample_page - 1) * nrow + 1
    sample_end <- min(sample_page * nrow, num_samples)
    NMRExperiment <- NMRExperiment[sample_start:sample_end]
    region_start <- (region_page - 1) * ncol + 1
    region_end <- min(region_page * ncol, num_regions)
    chemshift_range <- chemshift_range[region_start:region_end]

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
        ggplot2::scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_si("")))

    if (length(chemshift_range) > 1) {
        gplt <- gplt +
            ggplot2::facet_grid(
                rows = ggplot2::vars(.data$NMRExperiment),
                cols = ggplot2::vars(.data$region),
                scales = "free_x"
            )
    } else {
        gplt <- gplt +
            ggplot2::facet_grid(
                rows = ggplot2::vars(.data$NMRExperiment),
                scales = "free_x"
            )
    }
    gplt
}
