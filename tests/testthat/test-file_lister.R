test_that("file_lister lists files matching the glob, as plain character paths", {
    tmp_dir <- withr::local_tempdir()
    file.create(file.path(tmp_dir, c("Sample10.zip", "Sample20.zip", "notes.txt")))

    result <- file_lister(tmp_dir, "*.zip")

    expect_type(result, "character")
    expect_setequal(
        result,
        file.path(tmp_dir, c("Sample10.zip", "Sample20.zip"))
    )
})

test_that("file_lister returns an empty character vector when nothing matches", {
    tmp_dir <- withr::local_tempdir()
    file.create(file.path(tmp_dir, "notes.txt"))

    result <- file_lister(tmp_dir, "*.zip")

    expect_equal(result, character(0))
})

test_that("file_lister works on the package's demo dataset directory", {
    dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")

    result <- file_lister(dir_to_demo_dataset, "*.zip")

    expect_length(result, 3)
    expect_true(all(grepl("\\.zip$", result)))
})
