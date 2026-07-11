test_that("nmr_meta_get works", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = c(0:10),
        data_1r = matrix(sample(0:43, replace = TRUE), nrow = 4),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30", "40")))
    )
    dataset[["metadata"]][["external"]][["NMRExperiment"]] <- as.character(dataset[["metadata"]][["external"]][["NMRExperiment"]])
    meta <- nmr_meta_get(dataset, groups = "external")
    expect_equal(meta[[1, 1]], "10")
})

test_that("nmr_meta_get errors on an unknown group name", {
    dataset <- readRDS(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))

    # Real group names are present (sanity check the fixture is as expected).
    expect_true(all(c("external", "info") %in% nmr_meta_groups(dataset)))

    expect_error(
        nmr_meta_get(dataset, groups = "nonexistent_group_xyz"),
        regexp = "missing groups"
    )
    expect_error(
        nmr_meta_get(dataset, groups = "nonexistent_group_xyz"),
        regexp = "nonexistent_group_xyz"
    )
})

test_that("nmr_meta_get works for known groups", {
    dataset <- readRDS(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
    expect_no_error(meta <- nmr_meta_get(dataset, groups = "external"))
    expect_true("NMRExperiment" %in% colnames(meta))
})

test_that("nmr_meta_add_tidy_excel works", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = c(0:10),
        data_1r = matrix(sample(0:43, replace = TRUE), nrow = 4),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30", "40")))
    )
    dataset[["metadata"]][["external"]][["NMRExperiment"]] <- as.character(dataset[["metadata"]][["external"]][["NMRExperiment"]])

    dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
    MeOH_plasma_extraction_xlsx <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")

    dataset <- nmr_meta_add_tidy_excel(dataset, MeOH_plasma_extraction_xlsx)
    expect_match(dataset[["metadata"]][["external"]][["SubjectID"]][[1]], "Ana", ignore.case = TRUE)
})


test_that("nmr_meta_add works", {
    dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
    dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
    MeOH_plasma_extraction_xlsx <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
    exp_subj_id <- readxl::read_excel(MeOH_plasma_extraction_xlsx, sheet = 1)
    dataset <- nmr_meta_add(dataset, metadata = exp_subj_id, by = "NMRExperiment")
    meta <- nmr_meta_get(dataset, groups = "external")
    expect_match(meta[[1, 2]], "Ana", ignore.case = TRUE)
})

test_that("nmr_meta_add accepts a composite (multi-column) unnamed join key", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 0:10,
        data_1r = matrix(sample(0:43, 22, replace = TRUE), nrow = 2),
        metadata = list(external = data.frame(
            NMRExperiment = c("1", "2"),
            K1 = c("a", "b"),
            K2 = c("x", "y"),
            conflictCol = c("same1", "same2"),
            stringsAsFactors = FALSE
        ))
    )

    metadata_df <- data.frame(
        K1 = c("a", "b"),
        K2 = c("x", "y"),
        conflictCol = c("same1", "same2"), # identical values -> not a real conflict
        NewCol = c("n1", "n2"),
        stringsAsFactors = FALSE
    )

    expect_no_error(
        result <- nmr_meta_add(dataset, metadata_df, by = c("K1", "K2"))
    )

    meta <- nmr_meta_get(result, groups = "external")
    expect_true("NewCol" %in% colnames(meta))
    expect_false("K2__REMOVE__" %in% colnames(meta))
    expect_equal(meta$NewCol, c("n1", "n2"))
    # The (non-conflicting) shared column should be untouched.
    expect_equal(meta$conflictCol, c("same1", "same2"))
})

test_that("nmr_meta_add accepts a composite named join key (by = c(left = right))", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 0:10,
        data_1r = matrix(sample(0:43, 22, replace = TRUE), nrow = 2),
        metadata = list(external = data.frame(
            NMRExperiment = c("1", "2"),
            K1 = c("a", "b"),
            K2 = c("x", "y"),
            conflictCol = c("same1", "same2"),
            stringsAsFactors = FALSE
        ))
    )

    metadata_df <- data.frame(
        K1 = c("a", "b"),
        K2 = c("x", "y"),
        conflictCol = c("same1", "same2"),
        NewCol = c("n1", "n2"),
        stringsAsFactors = FALSE
    )

    expect_no_error(
        result <- nmr_meta_add(dataset, metadata_df, by = c(K1 = "K1", K2 = "K2"))
    )

    meta <- nmr_meta_get(result, groups = "external")
    expect_true("NewCol" %in% colnames(meta))
    expect_equal(meta$NewCol, c("n1", "n2"))
})

test_that("nmr_meta_add reports genuine column conflicts with a composite join key", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 0:10,
        data_1r = matrix(sample(0:43, 22, replace = TRUE), nrow = 2),
        metadata = list(external = data.frame(
            NMRExperiment = c("1", "2"),
            K1 = c("a", "b"),
            K2 = c("x", "y"),
            conflictCol = c("same1", "same2"),
            stringsAsFactors = FALSE
        ))
    )

    metadata_df <- data.frame(
        K1 = c("a", "b"),
        K2 = c("x", "y"),
        conflictCol = c("DIFFERENT1", "DIFFERENT2"), # genuinely conflicting values
        stringsAsFactors = FALSE
    )

    expect_error(
        nmr_meta_add(dataset, metadata_df, by = c("K1", "K2")),
        regexp = "conflictCol"
    )
})

test_that("nmr_meta_export warns with a comma-separated list of missing groups", {
    dataset <- readRDS(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
    xlsx_file <- withr::local_tempfile(fileext = ".xlsx")

    expect_warning(
        nmr_meta_export(
            dataset,
            xlsx_file,
            groups = c("nonexistent_group_one", "nonexistent_group_two", "external")
        ),
        regexp = "nonexistent_group_one, nonexistent_group_two"
    )

    # The export should still succeed with the remaining valid group.
    expect_true(file.exists(xlsx_file))
})

test_that("nmr_meta_export does not warn when all requested groups are present", {
    dataset <- readRDS(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
    xlsx_file <- withr::local_tempfile(fileext = ".xlsx")

    expect_no_warning(
        nmr_meta_export(dataset, xlsx_file, groups = "external")
    )
    expect_true(file.exists(xlsx_file))
})
