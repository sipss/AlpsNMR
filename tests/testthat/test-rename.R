test_that("test renaming of samples works", {
    # Create a dataset with two samples
    ds <- new_nmr_dataset_1D(
        ppm_axis = 1, 
        data_1r = matrix(1, nrow = 2), 
        metadata = list(external = data.frame(NMRExperiment = c("10", "20")))
    )
    # Rename sample 20 to newname:
    ds2 <- rename(ds, "newname" = "20")
    # Replace nmr_meta_get_column with names() once the behaviour change in names() happens
    expect_equal(nmr_meta_get_column(ds2), c("10", "newname"))
})
