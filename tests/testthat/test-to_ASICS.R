context("test-to_ASICS")

test_that("to_ASICS", {
dataset <- new_nmr_dataset_1D(ppm_axis = c(0:10),
                              data_1r = matrix(sample(0:43,replace = FALSE), nrow = 4),
                              metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30","40"))))
asics = to_ASICS(dataset)
expect_true(is.integer((asics@spectra@Dim)))
})

test_that("startsWith",{
startsWith <- function(x, prefix) {
  return(substring(x, 1, nchar(prefix)) == prefix)
}
expect_true(startsWith(22,22))
})