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
#' stopifnot(as.numeric(peak2peak_dist) == c(6, 0.1, 2, 0.9, 1, 6))
nmr_get_peak_distances <- function(peak_data, same_sample_dist_factor = 3) {
    peak_matrix <- matrix(peak_data$ppm, ncol = 1)
    rownames(peak_matrix) <- peak_data$peak_id
    peak2peak_dist <- peak2peak_distance(peak_matrix, distance_method = "euclidean")
    peak_groups <- peak_data |>
        dplyr::select(.data$NMRExperiment, .data$peak_id) |>
        dplyr::group_by(.data$NMRExperiment) %>%
        dplyr::summarise(peak_groups = list(.data$peak_id)) |>
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
#' @param num_clusters If you want to fix the number of clusters. Leave it to `NULL` to use a rough estimation
#'
#' @return A list including: The `peak_data` with an
#' additional "cluster" column, cluster (the hierarchical cluster),
#'  num_clusters (an estimation of the number of clusters)
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
nmr_peak_clustering <- function(peak_data, peak2peak_dist = NULL, num_clusters = NULL) {
    if (is.null(peak2peak_dist)) {
        peak2peak_dist <- nmr_get_peak_distances(peak_data)
    }
    cluster <- stats::hclust(d = peak2peak_dist, method = "complete")
    if (is.null(num_clusters)) {
        # FIXME: Implement something more robust to estimate num_clusters or height to cut:
        num_clusters <- which.max(-diff(sort(cluster$height, decreasing = TRUE)))+1
    }
    peak_data$cluster <- stats::cutree(cluster, k = num_clusters)
    list(
        peak_data = peak_data,
        cluster = cluster,
        num_clusters = num_clusters
    )
}


#' Build a peak table from the clustered peak list
#'
#' @param peak_data A peak list, with the cluster column
#'
#' @return A matrix, with one row per sample and one column per peak
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
nmr_build_peak_table <- function(peak_data) {
    if (!"cluster" %in% colnames(peak_data)) {
        stop("Please run nmr_peak_clustering() first")
    }
    peak_data <- peak_data |>
        dplyr::group_by(.data$cluster) |>
        dplyr::mutate(ppm_ref = round(stats::median(.data$ppm), 4)) |>
        dplyr::ungroup()
    
    peak_table <- peak_data |>
        dplyr::select(.data$NMRExperiment, .data$ppm_ref, .data$area) |>
        tidyr::pivot_wider(names_from = "ppm_ref", values_from = "area") |>
        tibble::column_to_rownames("NMRExperiment") |>
        as.matrix()
    
    nmr_exp <- stringr::str_sort(unique(peak_data$NMRExperiment), numeric = TRUE)
    
    peak_table[nmr_exp, , drop = FALSE]
}