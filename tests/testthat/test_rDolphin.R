test_that("to rDolphin files works", {
skip_on_bioc()
skip_if_not_installed("rDolphin")
dataset <- new_nmr_dataset_1D(ppm_axis = 1:10,
                              data_1r = matrix(sample(-1:88,replace = FALSE), nrow = 9),
                              metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90"))))
dataset[["metadata"]][["external"]][["SubjectID"]] = c("10", "20", "30", "40", "50", "60", "70", "80", "90")
dataset[["metadata"]][["external"]][["Group"]] = c("a", "a", "a", "a", "a", "b", "b", "b", "b")

blood = files_to_rDolphin(dataset, "blood")
cell = files_to_rDolphin(dataset, "cell")
urine = files_to_rDolphin(dataset, "urine")

expect_true(is.list(blood))
expect_true(is.list(cell))
expect_true(is.list(urine))

})
