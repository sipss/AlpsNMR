test_that("nmr_identify_regions_blood", {
  ppm_to_assign <- c(4.06)
  assignation <- nmr_identify_regions_blood(ppm_to_assign)
  expect_true(is.data.frame(assignation))
})

test_that("nmr_identify_regions_urine ", {
  ppm_to_assign <- c(4.06)
  assignation <- nmr_identify_regions_urine (ppm_to_assign)
  expect_true(is.data.frame(assignation))
})

test_that("nmr_identify_regions_cell  ", {
  ppm_to_assign <- c(4.06)
  assignation <- nmr_identify_regions_cell (ppm_to_assign)
  expect_true(is.data.frame(assignation))
})

