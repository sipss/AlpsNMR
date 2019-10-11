context("test-peak_detection")

test_that("nmr_detect_peaks & nmr_align_find_ref & nmr_align works", {
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
  dataset <- nmr_interpolate_1D(dataset, axis = c(min = 1, max = 2, by = 0.002))
  peak_table <- nmr_detect_peaks(dataset,
                                 nDivRange_ppm = 0.1,
                                 scales = seq(1, 16, 2),
                                 baselineThresh = 1e+03, SNR.Th = 2)
  NMRExp_ref <- nmr_align_find_ref(dataset, peak_table)
  dataset <- nmr_align(dataset, 
                       peak_table, 
                       NMRExp_ref, 
                       maxShift_ppm = 0.0015, 
                       acceptLostPeak = FALSE)
  expect_true(is.integer(dim(peak_table)))
  expect_true(is.character(NMRExp_ref))
  expect_true(is.matrix(dataset[["data_1r"]]))
  expect_true(is.numeric(dataset[["data_1r"]][[1]]))
})

