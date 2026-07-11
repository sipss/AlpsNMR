# Regression tests for bugs found and fixed in R/nmr_meta.R by a prior
# critical code review. Each test_that block below is written to fail
# against the pre-fix code and pass against the current (fixed) code.

# ---------------------------------------------------------------------------
# 1. nmr_meta_add(): multi-column join key (`by`) handling
#
# Previously:
#   by_left <- ifelse(is.null(names(by)), by, names(by))
# `ifelse()` returns a vector the same length as its *test*. Since
# `is.null(names(by))` is always a length-1 logical, `by_left` was silently
# truncated to a single element (the first join column) whenever `by` had
# more than one element - whether or not `by` was named. As a consequence,
# the second (and later) join key(s) were NOT excluded from
# `existing_vars <- setdiff(colnames(nmr_meta), by_left)`, so the join key
# itself was mistakenly treated as a candidate "conflicting" column. Because
# a genuine join key never gets a `__REMOVE__` suffix column from
# `dplyr::left_join()`, the subsequent `identical(nmr_meta_new[[col1]],
# nmr_meta_new[[col2]])` check compared the key column against `NULL`,
# which is never identical, and `nmr_meta_add()` aborted with a spurious
# "column conflict" error even though there was no real conflict.
#
# The fix uses a scalar `if`/`else` so `by_left` always contains *all* of
# the join key names, not just the first.
# ---------------------------------------------------------------------------

test_that("nmr_meta_add works with a composite (multi-column) unnamed join key", {
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

    # Before the fix, only "K1" was excluded from the conflict check, so "K2"
    # (a legitimate join key with no "__REMOVE__" counterpart) was wrongly
    # compared against NULL and this call raised
    # "Can't add metadata because of column conflict at: K2".
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

test_that("nmr_meta_add works with a composite named join key (by = c(left = right))", {
    # Same scenario as above but exercising the `names(by)` branch of the
    # `by_left` computation directly (the branch that held the original bug).
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

    # by = c(K1 = "K1", K2 = "K2") is a *named* multi-element vector, the
    # exact shape that triggered the truncation bug in
    # `ifelse(is.null(names(by)), by, names(by))` (names(by) got truncated
    # to names(by)[1] == "K1").
    expect_no_error(
        result <- nmr_meta_add(dataset, metadata_df, by = c(K1 = "K1", K2 = "K2"))
    )

    meta <- nmr_meta_get(result, groups = "external")
    expect_true("NewCol" %in% colnames(meta))
    expect_equal(meta$NewCol, c("n1", "n2"))
})

test_that("nmr_meta_add still reports genuine column conflicts with a composite join key", {
    # Sanity check that the fix didn't disable conflict detection altogether:
    # a real conflict (differing values in a shared, non-key column) must
    # still raise an error.
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

# ---------------------------------------------------------------------------
# 2. nmr_meta_get(): unknown `groups=` values must now error
#
# Previously an unrecognized group name in `groups=` was silently ignored:
# `metadata_list[groups]` with a name not present in `metadata_list` just
# produced an `NA`-named element that purrr::map/flatten quietly dropped,
# so the function returned only "NMRExperiment" with no error at all. The
# fix mirrors the existing `columns=` validation and now calls
# `rlang::abort()` listing the offending group name(s).
# ---------------------------------------------------------------------------

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

test_that("nmr_meta_get still works for known groups (no false positives)", {
    dataset <- readRDS(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
    expect_no_error(meta <- nmr_meta_get(dataset, groups = "external"))
    expect_true("NMRExperiment" %in% colnames(meta))
})

# ---------------------------------------------------------------------------
# 3. nmr_meta_export(): missing-groups warning message formatting
#
# Previously:
#   warning("...:\n", paste(groups[!groups_present]), collapse = ", ")
# `collapse` was passed as an (unused/invalid) argument to `warning()`
# itself, not to `paste()`, so `paste()` fell back to its default
# `sep = " "` with no `collapse`, and multiple missing group names were not
# comma-separated. The fix moves `collapse = ", "` into the `paste()` call.
# ---------------------------------------------------------------------------

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
