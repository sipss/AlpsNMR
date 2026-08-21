## peak2peak_distance ------------------------------------------------------

test_that("peak2peak_distance delegates to stats::dist for its documented methods", {
    x <- matrix(c(1, 10, 2, 20, 4, 50, 0, 0), ncol = 2, byrow = TRUE)
    rownames(x) <- paste0("Peak", 1:4)

    for (method in c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
        expect_equal(
            as.numeric(peak2peak_distance(x, distance_method = method)),
            as.numeric(stats::dist(x, method = method))
        )
    }
})

test_that("peak2peak_distance's sd_scaled_euclidean scales each column by its standard deviation first", {
    x <- matrix(c(1, 10, 2, 20, 4, 50, 0, 0), ncol = 2, byrow = TRUE)
    rownames(x) <- paste0("Peak", 1:4)

    d <- peak2peak_distance(x, distance_method = "sd_scaled_euclidean")
    expected <- stats::dist(scale(x, center = FALSE, scale = TRUE), method = "euclidean")
    expect_equal(as.numeric(d), as.numeric(expected))
})

test_that("peak2peak_distance rejects an unsupported distance method", {
    x <- matrix(1:4, ncol = 2)
    expect_error(
        peak2peak_distance(x, distance_method = "not-a-real-distance"),
        "Unsupported distance"
    )
})

test_that("mahalanobis_distance whitens by the sample covariance before computing Euclidean distances", {
    x <- matrix(c(1, 2, 2, 1, 4, 5, 0, 0), ncol = 2, byrow = TRUE)
    rownames(x) <- paste0("Peak", 1:4)

    d <- mahalanobis_distance(x)

    covinv <- solve(stats::cov(x))
    n <- nrow(x)
    expected <- matrix(0, n, n)
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            diff <- x[i, ] - x[j, ]
            expected[i, j] <- sqrt(as.numeric(t(diff) %*% covinv %*% diff))
        }
    }
    expect_equal(as.numeric(d), as.numeric(stats::as.dist(expected)))
})

## set_peak_distances_within_groups -----------------------------------------

test_that("set_peak_distances_within_groups overrides only within-group distances, and zeroes the diagonal", {
    m <- matrix(
        c(
            0, 1, 2, 3,
            1, 0, 4, 5,
            2, 4, 0, 6,
            3, 5, 6, 0
        ),
        nrow = 4, byrow = TRUE,
        dimnames = list(paste0("P", 1:4), paste0("P", 1:4))
    )
    d <- stats::as.dist(m)
    groups <- list(c("P1", "P2"), c("P3", "P4"))

    out <- as.matrix(set_peak_distances_within_groups(d, groups, value = 99))

    expect_equal(out["P1", "P2"], 99)
    expect_equal(out["P3", "P4"], 99)
    # Cross-group distances are untouched:
    expect_equal(out["P1", "P3"], 2)
    expect_equal(out["P1", "P4"], 3)
    expect_equal(out["P2", "P3"], 4)
    expect_equal(out["P2", "P4"], 5)
    # The diagonal stays zero:
    expect_equal(unname(diag(out)), rep(0, 4))
})

test_that("set_peak_distances_within_groups defaults to Inf", {
    m <- matrix(c(0, 1, 1, 0), nrow = 2, dimnames = list(c("P1", "P2"), c("P1", "P2")))
    out <- as.matrix(set_peak_distances_within_groups(stats::as.dist(m), list(c("P1", "P2"))))
    expect_equal(out["P1", "P2"], Inf)
})

test_that("set_peak_distances_within_groups masks the same pairs as a dense as.matrix() round-trip, without densifying", {
    # The implementation assigns directly into the triangular dist vector
    # (via index arithmetic) rather than round-tripping through as.matrix(),
    # to avoid doubling peak memory use; this checks the two give identical
    # results on a case too big to eyeball by hand.
    set.seed(42)
    n <- 40
    peak_names <- paste0("P", seq_len(n))
    d <- stats::dist(matrix(stats::runif(n), ncol = 1, dimnames = list(peak_names, NULL)))
    groups <- split(peak_names, sample(1:6, n, replace = TRUE))

    mask_dense_reference <- function(d, groups, value) {
        m <- as.matrix(d)
        for (ids in groups) {
            m[ids, ids] <- value
            diag(m) <- 0
        }
        stats::as.dist(m)
    }

    expected <- mask_dense_reference(d, groups, 999)
    actual <- set_peak_distances_within_groups(d, groups, value = 999)
    expect_equal(as.numeric(actual), as.numeric(expected))
    expect_equal(attr(actual, "Labels"), attr(d, "Labels"))
    expect_equal(attr(actual, "Size"), attr(d, "Size"))
    expect_s3_class(actual, "dist")
})

test_that("set_peak_distances_within_groups leaves singleton groups untouched", {
    m <- matrix(c(0, 1, 1, 0), nrow = 2, dimnames = list(c("P1", "P2"), c("P1", "P2")))
    out <- set_peak_distances_within_groups(stats::as.dist(m), list("P1", "P2"), value = 999)
    expect_equal(as.numeric(out), 1)
})

## split_peaks_into_ppm_components --------------------------------------------

test_that("split_peaks_into_ppm_components cuts wherever a gap exceeds the threshold", {
    peak_data <- data.frame(
        peak_id = paste0("Peak", 1:5),
        ppm = c(1.000, 1.001, 1.002, 2.000, 2.001)
    )
    # Gaps (ppb): 1, 1, 998, 1 -- only the third exceeds a 500 ppb threshold.
    components <- split_peaks_into_ppm_components(peak_data, max_dist_thresh_ppb = 500)

    expect_length(components, 2)
    expect_setequal(components[[1]]$peak_id, c("Peak1", "Peak2", "Peak3"))
    expect_setequal(components[[2]]$peak_id, c("Peak4", "Peak5"))
})

test_that("split_peaks_into_ppm_components keeps everything in one group when no gap exceeds the threshold", {
    peak_data <- data.frame(peak_id = paste0("Peak", 1:3), ppm = c(1, 1.001, 1.002))
    components <- split_peaks_into_ppm_components(peak_data, max_dist_thresh_ppb = 1e6)
    expect_length(components, 1)
    expect_equal(nrow(components[[1]]), 3)
})

test_that("split_peaks_into_ppm_components handles 0 and 1 peaks", {
    peak_data0 <- data.frame(peak_id = character(0), ppm = numeric(0))
    components0 <- split_peaks_into_ppm_components(peak_data0, max_dist_thresh_ppb = 100)
    expect_length(components0, 1)
    expect_equal(nrow(components0[[1]]), 0)

    peak_data1 <- data.frame(peak_id = "Peak1", ppm = 1.5)
    components1 <- split_peaks_into_ppm_components(peak_data1, max_dist_thresh_ppb = 100)
    expect_length(components1, 1)
    expect_equal(components1[[1]]$peak_id, "Peak1")
})

## cluster_one_ppm_component ---------------------------------------------------

test_that("cluster_one_ppm_component returns a single offset cluster for a single peak, with no hclust tree", {
    comp_peak_data <- data.frame(NMRExperiment = "10", peak_id = "Peak1", ppm = 1.5, gamma_ppb = 100)
    res <- cluster_one_ppm_component(comp_peak_data, max_dist_thresh_ppb = 300, cluster_offset = 5L)

    expect_equal(res$peak_data$cluster, 6L)
    expect_null(res$cluster)
    expect_equal(res$num_clusters, 1L)
    expect_null(res$num_cluster_estimation)
})

test_that("cluster_one_ppm_component offsets cluster ids and matches a directly-clustered reference", {
    comp_peak_data <- data.frame(
        NMRExperiment = c("10", "10", "20", "20"),
        peak_id = paste0("Peak", 1:4),
        ppm = c(1, 2, 1.1, 2.2),
        gamma_ppb = 100
    )
    reference <- nmr_peak_clustering(comp_peak_data, num_clusters = 2)
    res <- cluster_one_ppm_component(comp_peak_data, max_dist_thresh_ppb = 300, cluster_offset = 10L)

    expect_equal(res$num_clusters, 2L)
    ref_map <- setNames(reference$peak_data$cluster, reference$peak_data$peak_id)
    res_map <- setNames(res$peak_data$cluster - 10L, res$peak_data$peak_id)
    same_partition <- function(a, b) {
        a <- a[order(names(a))]
        b <- b[order(names(b))]
        identical(as.integer(factor(a)), as.integer(factor(b)))
    }
    expect_true(same_partition(ref_map, res_map))
    expect_true(all(res$peak_data$cluster > 10L))
})

## get_max_dist_ppb_for_num_clusters -------------------------------------------

test_that("get_max_dist_ppb_for_num_clusters handles a length-1 num_clusters (cutree drops to a plain vector)", {
    peak_list <- data.frame(peak_id = paste0("Peak", 1:3), ppm = c(1, 1.05, 1.1))
    d <- stats::dist(matrix(peak_list$ppm, dimnames = list(peak_list$peak_id, NULL)))
    cluster <- stats::hclust(d, method = "complete")

    # A single candidate k: stats::cutree() returns a bare named vector (no
    # dim) in this case rather than a matrix, which used to break the
    # subsequent matrix-style indexing.
    out <- get_max_dist_ppb_for_num_clusters(3, peak_list, cluster, max_dist_thresh_ppb = 1000)
    expect_equal(out$num_clusters, 3)
    expect_equal(out$max_distance_ppb, 0)
})

## nmr_get_peak_distances -----------------------------------------------------

test_that("nmr_get_peak_distances computes euclidean ppm distances and inflates within-sample pairs", {
    peak_data <- data.frame(
        NMRExperiment = c("10", "10", "20", "20"),
        peak_id = paste0("Peak", 1:4),
        ppm = c(1, 2, 1.1, 3)
    )
    peak2peak_dist <- nmr_get_peak_distances(peak_data)

    # stats::dist() order: (Peak2,Peak1), (Peak3,Peak1), (Peak4,Peak1),
    # (Peak3,Peak2), (Peak4,Peak2), (Peak4,Peak3). Peak1/Peak2 (NMRExperiment
    # "10") and Peak3/Peak4 (NMRExperiment "20") are same-sample pairs, so
    # their raw distance (1 and 6, respectively) is replaced by
    # same_sample_dist_factor (3, the default) times the overall max raw
    # distance (|1 - 3| = 2), i.e. 6.
    expect_equal(as.numeric(peak2peak_dist), c(6, 0.1, 2, 0.9, 1, 6), tolerance = 1e-8)
})

test_that("nmr_get_peak_distances's same_sample_dist_factor scales the within-sample override", {
    peak_data <- data.frame(
        NMRExperiment = c("10", "10", "20"),
        peak_id = paste0("Peak", 1:3),
        ppm = c(1, 2, 10)
    )
    # Max raw distance is |1 - 10| = 9; same-sample pair is (Peak1, Peak2).
    d1 <- nmr_get_peak_distances(peak_data, same_sample_dist_factor = 2)
    d2 <- nmr_get_peak_distances(peak_data, same_sample_dist_factor = 5)

    m1 <- as.matrix(d1)
    m2 <- as.matrix(d2)
    expect_equal(m1["Peak1", "Peak2"], 2 * 9)
    expect_equal(m2["Peak1", "Peak2"], 5 * 9)
    # Cross-sample distances are unaffected by the factor:
    expect_equal(m1["Peak1", "Peak3"], m2["Peak1", "Peak3"])
})

## nmr_peak_clustering ---------------------------------------------------------

test_that("nmr_peak_clustering assigns clusters, estimates their number, and computes ppm_ref as the per-cluster median", {
    peak_data <- data.frame(
        NMRExperiment = c("10", "10", "20", "20"),
        peak_id = paste0("Peak", 1:4),
        ppm = c(1, 2, 1.1, 2.2),
        gamma_ppb = 100
    )
    clustering_result <- nmr_peak_clustering(peak_data)

    expect_true("cluster" %in% colnames(clustering_result$peak_data))
    expect_equal(clustering_result$num_clusters, 2L)
    expect_setequal(clustering_result$peak_data$cluster, c(1L, 2L))
    # Peak1 (ppm=1) and Peak3 (ppm=1.1) form one cluster, Peak2 (ppm=2) and
    # Peak4 (ppm=2.2) the other; ppm_ref is the median ppm within the cluster.
    expect_equal(nrow(clustering_result$wrong_clusters), 0L)
    expect_null(clustering_result$excluded_peaks)
    ppm_ref_by_peak <- setNames(clustering_result$peak_data$ppm_ref, clustering_result$peak_data$peak_id)
    expect_equal(unname(ppm_ref_by_peak["Peak1"]), median(c(1, 1.1)))
    expect_equal(unname(ppm_ref_by_peak["Peak2"]), median(c(2, 2.2)))
})

test_that("nmr_peak_clustering honours an explicit num_clusters and skips the estimation step", {
    peak_data <- data.frame(
        NMRExperiment = c("10", "10", "20", "20"),
        peak_id = paste0("Peak", 1:4),
        ppm = c(1, 2, 1.1, 2.2),
        gamma_ppb = 100
    )
    clustering_result <- nmr_peak_clustering(peak_data, num_clusters = 2)

    expect_equal(clustering_result$num_clusters, 2)
    expect_null(clustering_result$num_cluster_estimation)
})

test_that("nmr_peak_clustering's verbose flag reports the estimated max distance threshold", {
    peak_data <- data.frame(
        NMRExperiment = c("10", "10", "20", "20"),
        peak_id = paste0("Peak", 1:4),
        ppm = c(1, 2, 1.1, 2.2),
        gamma_ppb = 100
    )
    # max_dist_thresh_ppb defaults to signif(3 * median(gamma_ppb), 2) = 300
    expect_message(
        nmr_peak_clustering(peak_data, verbose = TRUE),
        "300 ppbs"
    )
})

test_that("nmr_peak_clustering's default (region-split) path matches clustering everything at once", {
    # 4 widely-separated ppm pairs (consecutive gaps of ~0.9-6.8 ppm, all far
    # past the default max_dist_thresh_ppb of 300 ppb = 0.3 ppm), each pair
    # from two different samples. The default auto-estimate path splits by
    # ppm and clusters each ppm-disjoint group independently (see
    # split_peaks_into_ppm_components), landing each pair in its own
    # single-cluster group; this checks that gives the same overall
    # partition (up to cluster-id relabeling) as forcing everything through
    # one global hclust tree via a pre-supplied peak2peak_dist.
    peak_data <- data.frame(
        NMRExperiment = rep(c("10", "20"), 4),
        peak_id = paste0("Peak", 1:8),
        ppm = c(1, 1.1, 2, 2.2, 8, 8.1, 9, 9.2),
        gamma_ppb = 100
    )

    fast <- nmr_peak_clustering(peak_data, verbose = TRUE)
    forced_global <- nmr_peak_clustering(peak_data, peak2peak_dist = nmr_get_peak_distances(peak_data))

    expect_equal(fast$num_clusters, 4L)
    expect_type(fast$cluster, "list")
    expect_length(fast$cluster, 4L) # four ppm-disjoint groups, one cluster each

    same_partition <- function(a, b) {
        a <- a[order(names(a))]
        b <- b[order(names(b))]
        identical(as.integer(factor(a)), as.integer(factor(b)))
    }
    fast_map <- setNames(fast$peak_data$cluster, fast$peak_data$peak_id)
    global_map <- setNames(forced_global$peak_data$cluster, forced_global$peak_data$peak_id)
    expect_true(same_partition(fast_map, global_map))
})

test_that("nmr_peak_clustering warns and excludes peaks when two peaks from the same sample land in the same cluster", {
    peak_data <- data.frame(
        NMRExperiment = c("10", "10", "20", "20"),
        peak_id = paste0("Peak", 1:4),
        ppm = c(1, 2, 1.1, 2.2),
        gamma_ppb = 100
    )
    # Forcing a single cluster necessarily merges both of NMRExperiment
    # "10"'s peaks (and both of "20"'s) into it, which is an ambiguous
    # (more than one peak per sample per cluster) assignment.
    expect_warning(
        clustering_result <- nmr_peak_clustering(peak_data, num_clusters = 1),
        "Ambiguity detected"
    )

    expect_equal(nrow(clustering_result$peak_data), 0L)
    expect_equal(nrow(clustering_result$wrong_clusters), 2L)
    expect_equal(nrow(clustering_result$excluded_peaks), 4L)
    expect_setequal(clustering_result$excluded_peaks$peak_id, peak_data$peak_id)
})

## nmr_build_peak_table --------------------------------------------------------

test_that("nmr_build_peak_table pivots the clustered peak list into a peak table matrix", {
    peak_data <- data.frame(
        NMRExperiment = c("10", "10", "20", "20"),
        peak_id = paste0("Peak", 1:4),
        ppm = c(1, 2, 1.1, 2.1),
        area = c(10, 20, 12, 22)
    )
    clustering_result <- nmr_peak_clustering(peak_data, num_clusters = 2)
    peak_table <- nmr_build_peak_table(clustering_result$peak_data)

    expect_s3_class(peak_table, "nmr_dataset_peak_table")
    m <- nmr_data(peak_table)
    expect_equal(dim(m), c(2L, 2L))
    expect_setequal(rownames(m), c("10", "20"))
    expect_equal(unname(m["10", ]), c(10, 20))
    expect_equal(unname(m["20", ]), c(12, 22))
})

test_that("nmr_build_peak_table requires nmr_peak_clustering() to have been run first", {
    peak_data <- data.frame(
        NMRExperiment = c("10", "20"),
        peak_id = c("Peak1", "Peak2"),
        ppm = c(1, 2),
        area = c(10, 20)
    )
    expect_error(
        nmr_build_peak_table(peak_data),
        "nmr_peak_clustering"
    )
})

test_that("nmr_build_peak_table orders rows and pulls extra metadata from the dataset argument", {
    peak_data <- data.frame(
        NMRExperiment = c("10", "10", "20", "20"),
        peak_id = paste0("Peak", 1:4),
        ppm = c(1, 2, 1.1, 2.1),
        area = c(10, 20, 12, 22)
    )
    clustering_result <- nmr_peak_clustering(peak_data, num_clusters = 2)

    dataset <- new_nmr_dataset_1D(
        ppm_axis = seq(0.5, 2.5, length.out = 10),
        data_1r = matrix(stats::runif(20), nrow = 2),
        # Deliberately listed in the opposite order to the peak table's
        # natural row order, and carrying an extra metadata column:
        metadata = list(external = data.frame(NMRExperiment = c("20", "10"), Extra = c("x", "y")))
    )

    peak_table <- nmr_build_peak_table(clustering_result$peak_data, dataset = dataset)

    expect_equal(rownames(nmr_data(peak_table)), c("20", "10"))
    meta <- nmr_meta_get(peak_table)
    expect_equal(meta$Extra, c("x", "y"))
})

## nmr_peak_clustering_plot ----------------------------------------------------

test_that("nmr_peak_clustering_plot requires exactly two NMRExperiments", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = seq(1, 3, length.out = 10),
        data_1r = matrix(stats::runif(20), nrow = 2),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20")))
    )
    peak_list_clustered <- data.frame(
        NMRExperiment = c("10", "20"),
        ppm = c(1.5, 1.5),
        cluster = c(1, 1),
        area = c(10, 11),
        intensity_raw = c(5, 5.5)
    )

    expect_error(
        nmr_peak_clustering_plot(
            dataset, peak_list_clustered,
            NMRExperiments = c("10", "20", "30"),
            chemshift_range = c(1, 3)
        ),
        "2 and only 2"
    )
})

test_that("nmr_peak_clustering_plot builds a ggplot pairing matched-cluster peaks between two experiments", {
    skip_if_not_installed("ggplot2")

    ppm_axis <- seq(1, 3, length.out = 50)
    dataset <- new_nmr_dataset_1D(
        ppm_axis = ppm_axis,
        data_1r = matrix(stats::runif(2 * length(ppm_axis)), nrow = 2),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20")))
    )
    peak_list_clustered <- data.frame(
        NMRExperiment = c("10", "10", "20"),
        ppm = c(1.5, 2.5, 1.5),
        # cluster 1 has a peak on both experiments (a matched pair);
        # cluster 2 only has a peak on experiment "10".
        cluster = c(1, 2, 1),
        area = c(10, 20, 11),
        intensity_raw = c(5, 6, 5.5)
    )

    gplt <- nmr_peak_clustering_plot(
        dataset, peak_list_clustered,
        NMRExperiments = c("10", "20"),
        chemshift_range = c(1, 3)
    )

    expect_s3_class(gplt, "ggplot")
})
