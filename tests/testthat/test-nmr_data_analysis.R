## Prepare demo dataset
prepare_dataset <- function() {
    # 12 artificial samples created based on the 3 demo samples
    MeOH_plasma_extraction_dir <- system.file("dataset-demo", package = "AlpsNMR")
    MeOH_plasma_extraction_xlsx <- file.path(MeOH_plasma_extraction_dir, "dummy_metadata.xlsx")
    exp_subj_id <- readxl::read_excel(MeOH_plasma_extraction_xlsx, sheet = 1)

    zip_files <- fs::dir_ls(MeOH_plasma_extraction_dir, glob = "*.zip")

    dataset <- nmr_read_samples(sample_names = zip_files)
    dataset <- nmr_meta_add(dataset, metadata = exp_subj_id, by = "NMRExperiment")
    dataset <- nmr_interpolate_1D(dataset, axis = c(min = 3.7, max = 4.5, by = 2.3E-4))
    # nmr_baseline_removal() is deprecated (rate-limited cli::cli_warn(),
    # .frequency = "regularly"); suppress it here since it isn't what this
    # test is about, and its emission isn't reliable within a single session.
    dataset <- suppressWarnings(nmr_baseline_removal(dataset, lambda = 6, p = 0.01))
    dataset <- nmr_normalize(dataset, method = "area")

    metadata <- nmr_meta_get(dataset, groups = "external")
    metadata$Group <- c("A", "B", "B")
    # Artificially create a larger dataset
    larger_metadata <- rbind(metadata, metadata, metadata, metadata, metadata)

    larger_metadata$NMRExperiment <- as.character(
        seq(from = 10, by = 10, length.out = nrow(larger_metadata))
    )
    data_matrix <- nmr_data(dataset)
    dataset <- new_nmr_dataset_1D(
        ppm_axis = dataset$axis,
        data_1r = rbind(data_matrix, data_matrix, data_matrix, data_matrix, data_matrix),
        metadata = list(external = larger_metadata)
    )
    dataset
}

## Dataset can be used

test_that("nmr_data_analysis works", {
    dataset <- prepare_dataset()
    methodology <- plsda_auroc_vip_method(ncomp = 2)
    set.seed(123L)
    out <- nmr_data_analysis(
        dataset,
        y_column = "Group",
        identity_column = NULL,
        external_val = list(iterations = 1, test_size = 0.25),
        internal_val = list(iterations = 2, test_size = 0.25),
        data_analysis_method = methodology
    )
    expect_false(is.null(out))
})
