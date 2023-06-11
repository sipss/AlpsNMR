test_that("nmr_dataset_autophase works", {
    skip_if_not_installed("NMRphasing")
    dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
    dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
    nmr_dataset_autophase(dataset, method="NLS")
})