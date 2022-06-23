# Mahalanobis distance:
# https://stats.stackexchange.com/a/81710/62083
mahalanobis_distance <- function(x) {
    covmat <- stats::cov(x)
    dec <- chol(covmat)
    tmp <- forwardsolve(t(dec), t(x))
    colnames(tmp) <- rownames(x)
    stats::dist(t(tmp))
}

peak2peak_distance <- function(peak_matrix, distance_method = "euclidean") {
    STATS_METHODS <- c("euclidean", "maximum", "manhattan", "canberra",
                       "binary", "minkowski")
    if (distance_method %in% STATS_METHODS) {
        peak2peak_dist <- stats::dist(peak_matrix, method = distance_method)
    } else if (distance_method == "sd_scaled_euclidean") {
        peak_matrix_scaled <- scale(peak_matrix, center = FALSE, scale = TRUE)
        peak2peak_dist <- stats::dist(peak_matrix_scaled, method = "euclidean")
    } else if (distance_method == "mahalanobis") {
        peak2peak_dist <- mahalanobis_distance(peak_matrix)
    } else {
        stop(sprintf("Unsupported distance %s", distance_method))
    }
    peak2peak_dist
}

#' Override peak distances to infinity
#'
#' This function receives a distance matrix and a list of peak groups. Each group
#' consists of peaks that should not be grouped as the same peak (for instance because
#' they belong to the same sample). For each group, we set the distance between
#' all its peaks to infinity.
#'
#' @noRd
#'
#' @param dist_matrix A square matrix, where `dist_matrix[i,j]` is the distance
#'  from peak `i` to peak `j`. The matrix must have as row names and column names
#'  unique peak names.
#' @param peak_groups A list, where each element is a character vector with peak names
#' @param value `Inf` by default, but you could set to any other value
#'
#' @return An object of class "dist". See [stats::dist].
#'
set_peak_distances_within_groups <- function(dist_matrix, peak_groups, value = Inf) {
    # Set distances from pairs of peaks belonging to the same sample to Inf,
    # so they are never in the same cluster
    dist_matrix <- as.matrix(dist_matrix)
    for (peak_ids in peak_groups) {
        for (peak_i in peak_ids) {
            dist_matrix[peak_i, peak_ids] <- value
            dist_matrix[peak_ids, peak_i] <- value
            dist_matrix[peak_i, peak_i] <- 0
        }
    }
    stats::as.dist(dist_matrix)
}


#' Compute peak to peak distances
#'
#' @param peak_data A peak list
#' @param same_sample_dist_factor The distance between two peaks from the same
#' sample are set to this factor multiplied by the maximum of all the peak distances
#'
#' @return A dist object with the peak2peak distances
#' @export
#' @examples 
#' peak_data <- data.frame(
#'   NMRExperiment = c("10", "10", "20", "20"),
#'   peak_id = paste0("Peak", 1:4),
#'   ppm = c(1, 2, 1.1, 3)
#' )
#' peak2peak_dist <- nmr_get_peak_distances(peak_data)
#' stopifnot(abs(as.numeric(peak2peak_dist) - c(6, 0.1, 2, 0.9, 1, 6)) < 1E-8)
nmr_get_peak_distances <- function(peak_data, same_sample_dist_factor = 3) {
    peak_matrix <- matrix(peak_data$ppm, ncol = 1)
    rownames(peak_matrix) <- peak_data$peak_id
    peak2peak_dist <- peak2peak_distance(peak_matrix, distance_method = "euclidean")
    peak_groups <- peak_data %>%
        dplyr::select(.data$NMRExperiment, .data$peak_id) %>%
        dplyr::group_by(.data$NMRExperiment) %>%
        dplyr::summarise(peak_groups = list(.data$peak_id)) %>%
        tibble::deframe()
    max_dist <- max(peak2peak_dist)
    peak2peak_dist <- set_peak_distances_within_groups(
        dist_matrix = peak2peak_dist,
        peak_groups = peak_groups,
        value = same_sample_dist_factor*max_dist
    )
    peak2peak_dist
}

#' Peak clustering
#'
#' @param peak_data The peak list
#' @param peak2peak_dist The distances obtained with [nmr_get_peak_distances].
#'  If NULL it is computed from `peak_data`
#' @param num_clusters If you want to fix the number of clusters. Leave `NULL` if you want to estimate it
#' @param max_dist_thresh_ppb To estimate the number of clusters, we enforce a limit on how far two peaks of the same cluster may be. By default
#' this threshold will be computed as 3 times the median peak width (gamma), as given in the peak list.
#' @param verbose A logical vector to print additional information
#' @return A list including:
#'  - The `peak_data` with an additional "cluster" column
#'  - cluster: the hierarchical cluster
#'  -  num_clusters: an estimation of the number of clusters
#'  -  num_cluster_estimation: A list with tables and plots to justify the number of cluster estimation
#' @export
#'
#' @examples
#' peak_data <- data.frame(
#'   NMRExperiment = c("10", "10", "20", "20"),
#'   peak_id = paste0("Peak", 1:4),
#'   ppm = c(1, 2, 1.1, 3)
#' )
#' clustering_result <- nmr_peak_clustering(peak_data)
#' peak_data <- clustering_result$peak_data
#' stopifnot("cluster" %in% colnames(peak_data))
nmr_peak_clustering <- function(peak_data, peak2peak_dist = NULL, num_clusters = NULL, max_dist_thresh_ppb = NULL, verbose = FALSE) {
    if (is.null(peak2peak_dist)) {
        peak2peak_dist <- nmr_get_peak_distances(peak_data)
    }
    num_cluster_estimation <- NULL
    cluster <- stats::hclust(d = peak2peak_dist, method = "complete")
    if (is.null(num_clusters)) {
        if (is.null(max_dist_thresh_ppb)) {
            max_dist_thresh_ppb <- signif(3*stats::median(peak_data$gamma_ppb), digits = 2)
            if (verbose) {
                rlang::inform(c("i" = glue("The maximum distance between two peaks in the same cluster is of {max_dist_thresh_ppb} ppbs")))
            }
        }
        num_cluster_estimation <- estimate_num_clusters(
            peak_list = peak_data,
            cluster = cluster,
            max_dist_thresh_ppb = max_dist_thresh_ppb
        )
        num_clusters <- num_cluster_estimation$num_clusters
    }
    peak_data$cluster <- stats::cutree(cluster, k = num_clusters)
    
    # Estimate the ppm_ref
    # The digits are at least four. However if the max_dist_thresh is very small
    # then we need ppm references with more resolution as well.
    if (!is.null(max_dist_thresh_ppb)) {
        digits_for_ppmref <- min(4, 3 - floor(log10(max_dist_thresh_ppb)) + 1)
    } else {
        digits_for_ppmref <- 4
    }
    peak_data <- peak_data %>%
        dplyr::group_by(.data$cluster) %>%
        dplyr::mutate(ppm_ref = round(stats::median(.data$ppm), !!digits_for_ppmref)) %>%
        dplyr::ungroup()
    
    wrong_clusters <- peak_data %>%
        dplyr::group_by(.data$NMRExperiment, .data$ppm_ref) %>%
        dplyr::summarise(num_peaks = dplyr::n(), peak_ids = list(.data$peak_id), .groups = "drop") %>%
        dplyr::filter(.data$num_peaks > 1L)
    
    if (nrow(wrong_clusters) > 0) {
        wrong_peak_ids <- purrr::flatten_chr(wrong_clusters$peak_ids)
        rlang::warn(
            message = c(
                glue("Ambiguity detected in the peak clustering affecting {length(wrong_peak_ids)} out of {nrow(peak_data)} peaks in the dataset"),
                "i" = "Some samples have more than one peak in the same cluster.",
                "i" = "This may indicate either a too small number of clusters or some false peaks in the peak table",
                "i" = "As a workaround, we have removed all those problematic peaks",
                "i" = "You can recover some of the problematic peaks exploring the wrong_clusters or excluded_peaks, or decide to proceed without them."
            )
        )
        excluded_peaks <- dplyr::filter(peak_data, .data$peak_id %in% wrong_peak_ids)
        peak_data <- dplyr::filter(peak_data, !.data$peak_id %in% wrong_peak_ids)
    } else {
        excluded_peaks <- NULL
    }
    
    
    list(
        peak_data = peak_data,
        cluster = cluster,
        num_clusters = num_clusters,
        num_cluster_estimation = num_cluster_estimation,
        wrong_clusters = wrong_clusters,
        excluded_peaks = excluded_peaks
    )
}

get_max_dist_ppb_for_num_clusters <- function(num_clusters, peak_list, cluster) {
    peak_list$cluster <- stats::cutree(cluster, k = num_clusters)
    max_dist_ppm <- peak_list |>
        dplyr::group_by(.data$cluster) |> 
        dplyr::summarize(max_dist = max(.data$ppm) - min(.data$ppm), .groups = "drop") |>
        dplyr::pull("max_dist") |>
        max()
    1000*max_dist_ppm
}

estimate_num_clusters <- function(peak_list, cluster, max_dist_thresh_ppb) {
    peaks_per_sample <- peak_list |> 
        dplyr::group_by(.data$NMRExperiment) |> 
        dplyr::summarize(n = dplyr::n()) |>
        dplyr::pull("n")
    num_clusters_coarse <- seq.int(from = max(peaks_per_sample), to = sum(peaks_per_sample), by = 10)
    max_dist_coarse <- purrr::map_dbl(
        num_clusters_coarse,
        get_max_dist_ppb_for_num_clusters,
        peak_list = peak_list,
        cluster = cluster
    )
    num_clusters <- num_clusters_coarse[max_dist_coarse < max_dist_thresh_ppb][1]
    # Refine:
    num_clusters_fine <- c(
        seq.int(from = num_clusters - 19, to = num_clusters - 11),
        seq.int(from = num_clusters - 9, to = num_clusters - 1),
        seq.int(from = num_clusters + 1, to = num_clusters + 9)
    )
    max_dist_fine <- purrr::map_dbl(
        num_clusters_fine,
        get_max_dist_ppb_for_num_clusters,
        peak_list = peak_list,
        cluster = cluster
    )
    # Combine:
    num_clusters_vs_max_distance <- data.frame(
        num_clusters = c(num_clusters_coarse, num_clusters_fine),
        max_distance_ppb = c(max_dist_coarse, max_dist_fine)
    ) |> dplyr::arrange(num_clusters)
    num_clusters <- num_clusters_vs_max_distance |>
        dplyr::filter(.data$max_distance_ppb < !!max_dist_thresh_ppb) |>
        dplyr::pull("num_clusters")
    num_clusters <- num_clusters[1]
    gplt <- ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = .data$num_clusters, y = .data$max_distance_ppb), data = num_clusters_vs_max_distance, na.rm = TRUE) +
        ggplot2::geom_hline(yintercept = max_dist_thresh_ppb, color = "gray") + 
        ggplot2::geom_vline(xintercept = num_clusters, color = "red") +
        ggplot2::labs(x = "Number of clusters", y = "Max distance within cluster (ppb)")
    list(
        num_clusters = num_clusters,
        table = num_clusters_vs_max_distance,
        max_dist_thresh_ppb = max_dist_thresh_ppb,
        plot = gplt
    )
}

#' Plot clustering results
#'
#' @param dataset The [nmr_dataset_1D] object
#' @param peak_list_clustered A peak list table with a clustered column
#' @param NMRExperiments Two and only two experiments to compare in the plot
#' @param chemshift_range A region, make it so it does not cover a huge range (maybe 1ppm or less)
#'
#' @return A plot of the two experiments in the given chemshift range, with
#'  lines connecting peaks identified as the same and dots showing peaks without pairs
#' @export
#'
nmr_peak_clustering_plot <- function(dataset, peak_list_clustered, NMRExperiments, chemshift_range) {
    if (length(NMRExperiments) != 2) {
        rlang::abort("Please provide 2 and only 2 NMRExperiments")
    }
    
    spectra <- tidy(dataset, chemshift_range = chemshift_range, NMRExperiment = NMRExperiments)
    offset_for_plotting <- 0.2*diff(range(spectra$intensity))
    spec_rows_2 <- spectra$NMRExperiment == NMRExperiments[2]
    spectra$intensity[spec_rows_2] <- spectra$intensity[spec_rows_2] + offset_for_plotting
    peak_list_clustered2 <- dplyr::filter(
        peak_list_clustered,
        .data$NMRExperiment %in% !!NMRExperiments,
        .data$ppm > min(!!chemshift_range),
        .data$ppm < max(!!chemshift_range)
    )
    
    plc_rows_2 <- peak_list_clustered2$NMRExperiment == NMRExperiments[2]
    peak_list_clustered2$intensity[plc_rows_2] <- peak_list_clustered2$intensity[plc_rows_2] + offset_for_plotting
    
    for_segments <- dplyr::full_join(
        peak_list_clustered2 |>
            dplyr::filter(.data$NMRExperiment == !!NMRExperiments[1]) |>
            dplyr::select(.data$NMRExperiment, .data$ppm, .data$intensity, .data$cluster, .data$area),
        peak_list_clustered2 |>
            dplyr::filter(.data$NMRExperiment == !!NMRExperiments[2]) |>
            dplyr::select(.data$NMRExperiment, .data$ppm, .data$intensity, .data$cluster, .data$area),
        by = "cluster"
    )
    
    only_on_1 <- for_segments |> filter(is.na(.data$NMRExperiment.y))
    only_on_2 <- for_segments |> filter(is.na(.data$NMRExperiment.x))
    for_segments_12 <- for_segments |> filter(!is.na(.data$NMRExperiment.x), !is.na(.data$NMRExperiment.y))
    
    geom_txt <- get_geom_text()
    ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = .data$chemshift, y = .data$intensity, color = .data$NMRExperiment), data = spectra) +
        ggplot2::geom_hline(yintercept = c(0, offset_for_plotting), color = "gray") +

        ggplot2::geom_point(ggplot2::aes(x = .data$ppm.x, y = .data$intensity.x), data = only_on_1) + 
        geom_txt(ggplot2::aes(x = .data$ppm.x, y = .data$intensity.x, color = .data$NMRExperiment.x, label = signif(.data$area.x, 4)), data = only_on_1) + 

        ggplot2::geom_point(ggplot2::aes(x = .data$ppm.y, y = .data$intensity.y), data = only_on_2) +
        geom_txt(ggplot2::aes(x = .data$ppm.y, y = .data$intensity.y, color = .data$NMRExperiment.y, label = signif(.data$area.y, 4)), data = only_on_2) + 

        ggplot2::geom_segment(ggplot2::aes(x = .data$ppm.x, y = .data$intensity.x, xend = .data$ppm.y, yend = .data$intensity.y), data = for_segments_12) +
        geom_txt(ggplot2::aes(x = .data$ppm.x, y = .data$intensity.x, color = .data$NMRExperiment.x, label = signif(.data$area.x, 4)), data = for_segments_12) + 
        geom_txt(ggplot2::aes(x = .data$ppm.y, y = .data$intensity.y, color = .data$NMRExperiment.y, label = signif(.data$area.y, 4)), data = for_segments_12) + 
        ggplot2::scale_x_reverse() +
        ggplot2::scale_y_continuous(labels = scales::label_number (scale_cut = scales::cut_si(""))) +
        ggplot2::labs(x = "Chemical shift (ppm)", y = "Intensity")
}



#' Build a peak table from the clustered peak list
#'
#' @param peak_data A peak list, with the cluster column
#' @param dataset A [nmr_dataset_1D] object, to get the matadata
#' @return An [nmr_dataset_peak_table], containing the peak table and the annotations
#' @export
#' @examples
#' peak_data <- data.frame(
#'   NMRExperiment = c("10", "10", "20", "20"),
#'   peak_id = paste0("Peak", 1:4),
#'   ppm = c(1, 2, 1.1, 2.1),
#'   area = c(10, 20, 12, 22)
#' )
#' clustering_result <- nmr_peak_clustering(peak_data, num_clusters = 2)
#' peak_data <- clustering_result$peak_data
#' peak_table <- nmr_build_peak_table(peak_data)
#' stopifnot(ncol(peak_table) == 2)
nmr_build_peak_table <- function(peak_data, dataset) {
    if (!"cluster" %in% colnames(peak_data)) {
        stop("Please run nmr_peak_clustering() first")
    }
    
    peak_table <- peak_data %>%
        dplyr::select(.data$NMRExperiment, .data$ppm_ref, .data$area) %>%
        dplyr::mutate(ppm_ref = format(.data$ppm_ref)) %>%
        tidyr::pivot_wider(names_from = "ppm_ref", values_from = "area") %>%
        tibble::column_to_rownames("NMRExperiment") %>%
        as.matrix()
    
    nmr_exp <- stringr::str_sort(unique(peak_data$NMRExperiment), numeric = TRUE)
    new_nmr_dataset_peak_table(
        peak_table = peak_table[nmr_exp, , drop = FALSE],
        metadata = list(external = nmr_meta_get(dataset, groups = "external"))
    )
}