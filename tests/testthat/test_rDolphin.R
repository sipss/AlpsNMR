test_that("to rDolphin files works", {
    skip_on_bioc()
    skip_if_not_installed("rDolphin")
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:10,
        data_1r = matrix(sample(-1:88, replace = FALSE), nrow = 9),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90")))
    )
    dataset[["metadata"]][["external"]][["SubjectID"]] <- c("10", "20", "30", "40", "50", "60", "70", "80", "90")
    dataset[["metadata"]][["external"]][["Group"]] <- c("a", "a", "a", "a", "a", "b", "b", "b", "b")

    blood <- files_to_rDolphin(dataset, "blood")
    cell <- files_to_rDolphin(dataset, "cell")
    urine <- files_to_rDolphin(dataset, "urine")

    expect_true(is.list(blood))
    expect_true(is.list(cell))
    expect_true(is.list(urine))
})

test_that("files_to_rDolphin() reports the Group label to numeric code mapping", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:10,
        data_1r = matrix(
            stats::runif(9 * 10),
            nrow = 9
        ),
        metadata = list(
            external = data.frame(
                NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90"),
                stringsAsFactors = FALSE
            )
        )
    )
    dataset[["metadata"]][["external"]][["SubjectID"]] <-
        c("10", "20", "30", "40", "50", "60", "70", "80", "90")
    dataset[["metadata"]][["external"]][["Group"]] <-
        c("Control", "Control", "Control", "Control", "Control", "Patient", "Patient", "Patient", "Patient")

    messages <- testthat::capture_messages(
        result <- files_to_rDolphin(dataset, "blood")
    )
    mapping_message <- paste(messages, collapse = "\n")

    # The message must mention both original labels together with their
    # assigned numeric codes, so the user can see how the encoding was done.
    expect_match(mapping_message, "Control\\s*=\\s*1")
    expect_match(mapping_message, "Patient\\s*=\\s*2")

    # Sanity check that the returned data actually uses that same mapping.
    expect_equal(
        result[["meta_rDolphin"]][["type"]],
        c(1, 1, 1, 1, 1, 2, 2, 2, 2)
    )
})
