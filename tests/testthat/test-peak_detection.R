test_that("nmr_detect_peaks & nmr_align_find_ref & nmr_align & nmr_integrate_peak_position works", {
  dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
  dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
  data2 <- dataset[1:3]
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

  plot <- nmr_detect_peaks_plot(dataset,peak_table,NMRExp_ref)
  
  
  peak_table_integration = nmr_integrate_peak_positions(
  samples = dataset,
  peak_pos_ppm = list(c(2,3,4)),
  peak_width_ppm = NULL)
  expect_true(is.list(data2[["data_1r"]]))
  expect_true(is.numeric(peak_table_integration[["peak_table"]][[1]]))
  expect_true(is.integer(dim(peak_table)))
  expect_true(is.character(NMRExp_ref))
  expect_true(is.list(plot[1]))
  expect_true(is.matrix(dataset[["data_1r"]]))
  expect_true(is.numeric(dataset[["data_1r"]][[1]]))
})

test_that("nmr_integrate_regions works", {
  dataset <- new_nmr_dataset_1D(ppm_axis = 1:10,
                                data_1r = matrix(sample(0:99,replace = TRUE), nrow = 10),
                                metadata = list(external = data.frame(NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100"))))
  peak_table_integration = nmr_integrate_regions(
    samples = dataset,
    regions = list(random = c(2,5)),
    fix_baseline = TRUE)
  expect_true(is.matrix(peak_table_integration[["peak_table"]]))
  expect_true(is.numeric(peak_table_integration[["peak_table"]][[1]]))
})

