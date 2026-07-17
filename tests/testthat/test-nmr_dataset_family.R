build_family_test_dataset <- function() {
    metadata <- list(
        external = data.frame(NMRExperiment = c("10", "20", "30")),
        info = data.frame(NMRExperiment = c("10", "20", "30"), info_x = 1:3)
    )
    new_nmr_dataset(
        metadata,
        list(data_1r = list(runif(4), runif(4), runif(4))),
        list(list(1:4), list(1:4), list(1:4))
    )
}

## names.nmr_dataset_family ---------------------------------------------------

test_that("names.nmr_dataset_family returns the NMRExperiment column", {
    dataset <- build_family_test_dataset()
    expect_equal(names(dataset), c("10", "20", "30"))
})

## names<-.nmr_dataset_family --------------------------------------------------

test_that("names<- renames NMRExperiment consistently across every metadata table", {
    # Regression test for https://github.com/sipss/AlpsNMR/issues/62: the
    # reported workaround was to manually loop over every metadata table and
    # overwrite its NMRExperiment column. names<- must do exactly that in a
    # single call, so every table stays consistent with the new names.
    dataset <- build_family_test_dataset()

    names(dataset) <- c("A", "B", "C")

    expect_equal(names(dataset), c("A", "B", "C"))
    expect_equal(dataset$metadata$external$NMRExperiment, c("A", "B", "C"))
    expect_equal(dataset$metadata$info$NMRExperiment, c("A", "B", "C"))
})

test_that("names<- works on a single-sample dataset obtained by subsetting", {
    # Regression test for https://github.com/sipss/AlpsNMR/issues/62: renaming
    # a subsetted (single-sample) dataset used to fail with "number of items
    # to replace is not a multiple of replacement length".
    dataset <- build_family_test_dataset()
    single <- dataset[1]

    expect_no_warning(names(single) <- "MyNewName")
    expect_no_error(names(single) <- "MyNewName")

    expect_equal(names(single), "MyNewName")
    expect_equal(single$metadata$external$NMRExperiment, "MyNewName")
    expect_equal(single$metadata$info$NMRExperiment, "MyNewName")
})

test_that("names<- requires exactly one name per sample", {
    dataset <- build_family_test_dataset()
    expect_error(
        names(dataset) <- c("A", "B"),
        "length 3.*length 2"
    )
})

test_that("names<- rejects duplicated names", {
    dataset <- build_family_test_dataset()
    expect_error(
        names(dataset) <- c("A", "A", "B"),
        "should not be repeated"
    )
})
