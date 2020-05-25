# context("test-MUVR")
# 
# test_that("rdCV_PLS_RF & permutation_test_model work", {
#   # Too expensive test to run on CRAN/Bioconductor
#   skip_on_cran()
#   skip_on_bioc()
#   skip_if_not_installed("MUVR")
#   dataset <- new_nmr_dataset_1D(
#     ppm_axis = 1:10,
#     data_1r = matrix(sample(0:99,replace = TRUE), nrow = 10),
#     metadata = list(external = data.frame(
#       NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100"),
#       group = c("A", "A", "A", "A", "A", "B", "B", "B", "B", "B"), 
#       stringsAsFactors = FALSE))
#   )
#   meta <- nmr_meta_get(dataset, groups = "external")
#   model <- rdCV_PLS_RF(nmr_data(dataset),Y = meta$group, nOuter = 3, nInner = 2, nRep = 2, parallel = FALSE)
#   permutations = permutation_test_model(model, nPerm = 2, parallel = FALSE)
#   VIPs= model_VIP(model)
#   expect_true(is.numeric(model[["calcMins"]]))
#   expect_true(is.matrix(permutations))
#   expect_true(is.data.frame(VIPs))
# })

