test_that("create_sample_names strips zip inner paths before deduplicating", {
    sample_names <- c("foo.zip!/inner/sample1", "bar.zip!/inner/sample2")
    expect_equal(create_sample_names(sample_names), c("foo", "bar"))
})

test_that("create_sample_names falls back to vctrs::vec_as_names when still duplicated", {
    sample_names <- c("archive.zip!/inner/sample1", "archive.zip!/inner/sample2")
    expect_message(
        result <- create_sample_names(sample_names),
        "New names"
    )
    expect_false(anyDuplicated(result) > 0)
    expect_length(result, 2)
})

test_that("new_nmr_dataset builds a valid nmr_dataset object", {
    metadata <- list(external = data.frame(NMRExperiment = c("10", "20")))
    data_fields <- list(data_1r = list(runif(16), runif(32)))
    axis <- list(list(1:16), list(1:32))
    dataset <- new_nmr_dataset(metadata, data_fields, axis)

    expect_s3_class(dataset, "nmr_dataset")
    expect_s3_class(dataset, "nmr_dataset_family")
    expect_equal(dataset$num_samples, 2)
    expect_equal(dataset$metadata, metadata)
    expect_length(dataset$data_1r, 2)
    expect_length(dataset$data_1r[[1]], 16)
    expect_length(dataset$data_1r[[2]], 32)
    expect_equal(dataset$axis, axis)
})

test_that("new_nmr_dataset supports several data_ fields at once", {
    metadata <- list(external = data.frame(NMRExperiment = c("10", "20")))
    data_fields <- list(
        data_1r = list(runif(8), runif(8)),
        data_1i = list(runif(8), runif(8))
    )
    axis <- list(list(1:8), list(1:8))
    dataset <- new_nmr_dataset(metadata, data_fields, axis)

    expect_true(all(c("data_1r", "data_1i") %in% names(unclass(dataset))))
    expect_length(dataset$data_1i, 2)
})

test_that("is.nmr_dataset distinguishes nmr_dataset objects from other objects", {
    metadata <- list(external = data.frame(NMRExperiment = "10"))
    dataset <- new_nmr_dataset(metadata, list(data_1r = list(runif(4))), list(list(1:4)))

    expect_true(is.nmr_dataset(dataset))
    expect_false(is.nmr_dataset(list()))
    expect_false(is.nmr_dataset(1L))
})

test_that("[.nmr_dataset subsets metadata, data fields, axis and num_samples together", {
    metadata <- list(
        external = data.frame(NMRExperiment = c("10", "20", "30")),
        extra = data.frame(NMRExperiment = c("10", "20", "30"), value = c(1, 2, 3))
    )
    data_fields <- list(data_1r = list(runif(4), runif(5), runif(6)))
    axis <- list(list(1:4), list(1:5), list(1:6))
    dataset <- new_nmr_dataset(metadata, data_fields, axis)

    subset_dataset <- dataset[c(1, 3)]

    expect_equal(subset_dataset$num_samples, 2)
    expect_equal(subset_dataset$metadata$external$NMRExperiment, c("10", "30"))
    expect_equal(subset_dataset$metadata$extra$value, c(1, 3))
    expect_length(subset_dataset$data_1r, 2)
    expect_length(subset_dataset$data_1r[[1]], 4)
    expect_length(subset_dataset$data_1r[[2]], 6)
    expect_length(subset_dataset$axis, 2)
    expect_true(is.nmr_dataset(subset_dataset))
})

test_that("print.nmr_dataset and format.nmr_dataset describe the number of samples", {
    metadata <- list(external = data.frame(NMRExperiment = c("10", "20")))
    dataset <- new_nmr_dataset(metadata, list(data_1r = list(runif(2), runif(2))), list(list(1:2), list(1:2)))

    expect_equal(format(dataset), "An nmr_dataset (2 samples)")
    expect_output(print(dataset), "An nmr_dataset \\(2 samples\\)")
    expect_invisible(print(dataset))
})

test_that("validate_nmr_dataset accepts a valid object and returns it invisibly unchanged", {
    metadata <- list(external = data.frame(NMRExperiment = "10"))
    dataset <- new_nmr_dataset(metadata, list(data_1r = list(runif(4))), list(list(1:4)))
    expect_equal(validate_nmr_dataset(dataset), dataset)
})

test_that("validate_nmr_dataset rejects an object without the nmr_dataset class", {
    metadata <- list(external = data.frame(NMRExperiment = c("10", "20")))
    not_a_dataset <- structure(
        list(metadata = metadata, num_samples = 2),
        class = "nmr_dataset_family"
    )
    expect_error(validate_nmr_dataset(not_a_dataset), "Not an nmr_dataset object")
})

test_that("validate_nmr_dataset propagates nmr_dataset_family validation errors", {
    metadata <- list(external = data.frame(x = c("10", "20")))
    bad_dataset <- structure(
        list(metadata = metadata, num_samples = 2),
        class = c("nmr_dataset", "nmr_dataset_family")
    )
    expect_error(validate_nmr_dataset(bad_dataset), "NMRExperiment column")
})
