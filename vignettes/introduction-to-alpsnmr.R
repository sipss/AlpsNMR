## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  fig.width = 7,
  fig.height = 5,
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(AlpsNMR)
library(ggplot2)

## ------------------------------------------------------------------------
#plan(sequential()) # disable parallellization (default)
plan(multiprocess(workers = 4)) # enable parallellization with 4 workers

## ------------------------------------------------------------------------
MeOH_serum_extraction_dir <- system.file("dataset-demo", package = "AlpsNMR")
print(MeOH_serum_extraction_dir)

## ------------------------------------------------------------------------
fs::dir_ls(MeOH_serum_extraction_dir)

## ------------------------------------------------------------------------
MeOH_serum_extraction_xlsx <- file.path(MeOH_serum_extraction_dir, "dummy_metadata.xlsx")
exp_subj_id <- readxl::read_excel(MeOH_serum_extraction_xlsx, sheet = 1)
subj_id_age <- readxl::read_excel(MeOH_serum_extraction_xlsx, sheet = 2)
print(exp_subj_id)
print(subj_id_age)


## ----load-samples--------------------------------------------------------
zip_files <- fs::dir_ls(MeOH_serum_extraction_dir, glob = "*.zip")
print(zip_files)
dataset <- nmr_read_samples(sample_names = zip_files)
print(dataset)

## ------------------------------------------------------------------------
nmr_meta_get(dataset, groups = "external")

## ------------------------------------------------------------------------
nmr_meta_get(dataset, groups = "external")

## ------------------------------------------------------------------------
dataset <- nmr_meta_add(dataset, metadata = exp_subj_id, by = "NMRExperiment")
nmr_meta_get(dataset, groups = "external")

## ------------------------------------------------------------------------
dataset <- nmr_meta_add(dataset, metadata = subj_id_age, by = "SubjectID")
nmr_meta_get(dataset, groups = "external")

## ------------------------------------------------------------------------
ppm_res <- nmr_ppm_resolution(dataset)[[1]]
message("The ppm resolution is: ", format(ppm_res, digits = 2), " ppm")

## ------------------------------------------------------------------------
dataset <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))

## ------------------------------------------------------------------------
plot(dataset, NMRExperiment = c("10", "30"), chemshift_range = c(2.2, 2.8))

## ------------------------------------------------------------------------
regions_to_exclude <- list(water = c(4.6, 5), methanol = c(3.33, 3.39))
dataset <- nmr_exclude_region(dataset, exclude = regions_to_exclude)
plot(dataset, chemshift_range = c(4.2, 5.5))

## ------------------------------------------------------------------------
samples_10_20 <- filter(dataset, SubjectID == "Ana")
nmr_meta_get(samples_10_20, groups = "external")

## ------------------------------------------------------------------------
pca_outliers_rob <- nmr_pca_outliers_robust(dataset, ncomp = 3)
nmr_pca_outliers_plot(dataset, pca_outliers_rob)

## ------------------------------------------------------------------------
plot(dataset, chemshift_range = c(3.5,3.8))

## ------------------------------------------------------------------------
dataset = nmr_baseline_removal(dataset, lambda = 6, p = 0.01)
plot(dataset, chemshift_range = c(3.5,3.8))

## ------------------------------------------------------------------------
peak_table <- nmr_detect_peaks(dataset,
                               nDivRange_ppm = 0.1,
                               scales = seq(1, 16, 2),
                               baselineThresh = NULL, SNR.Th = 3)
NMRExp_ref <- nmr_align_find_ref(dataset, peak_table)
message("Your reference is NMRExperiment ", NMRExp_ref)
nmr_detect_peaks_plot(dataset, peak_table, NMRExperiment = "20", chemshift_range = c(3.5,3.8))

## ------------------------------------------------------------------------
nmr_exp_ref <- nmr_align_find_ref(dataset, peak_table)
dataset_align <- nmr_align(dataset, peak_table, nmr_exp_ref, maxShift_ppm = 0.0015, acceptLostPeak = FALSE)

## ------------------------------------------------------------------------
plot(dataset, chemshift_range = c(3.025, 3.063))
plot(dataset_align, chemshift_range = c(3.025, 3.063))

## ------------------------------------------------------------------------



## ------------------------------------------------------------------------
dataset_norm <- nmr_normalize(dataset_align, method = "pqn")

## ------------------------------------------------------------------------
diagnostic <- nmr_diagnose(dataset_norm)
print(diagnostic$norm_factor)
print(diagnostic$plot)

## ------------------------------------------------------------------------
peak_table_integration = nmr_integrate_peak_positions(
  samples = dataset_norm,
  peak_pos_ppm = peak_table$ppm,
  peak_width_ppm = 0.006)

peak_table_integration = get_integration_with_metadata(peak_table_integration)

## ------------------------------------------------------------------------
nmr_data(
  nmr_integrate_peak_positions(samples = dataset_norm,
                            peak_pos_ppm = c(4.1925, 4.183, 4.1775, 4.17),
                            peak_width_ppm = 0.006)
)

## ------------------------------------------------------------------------
pyroglutamic_acid_region <- c(4.15, 4.20)
plot(dataset_norm, chemshift_range = pyroglutamic_acid_region) +
  ggplot2::ggtitle("Pyroglutamic acid region")

## ------------------------------------------------------------------------
pyroglutamic_acid <- list(pyroglutamic_acid1 = c(4.19, 4.195),
                          pyroglutamic_acid2 = c(4.18, 4.186),
                          pyroglutamic_acid3 = c(4.175, 4.18),
                          pyroglutamic_acid4 = c(4.165, 4.172))
regions_basel_corr_ds <- nmr_integrate_regions(dataset_norm, pyroglutamic_acid, fix_baseline = TRUE)
regions_basel_corr_matrix <- nmr_data(regions_basel_corr_ds)
print(regions_basel_corr_matrix)


regions_basel_not_corr_ds <- nmr_integrate_regions(dataset_norm, pyroglutamic_acid, fix_baseline = FALSE)
regions_basel_not_corr_matrix <- nmr_data(regions_basel_not_corr_ds)
print(regions_basel_not_corr_matrix)

## ------------------------------------------------------------------------
dplyr::bind_rows(
  regions_basel_corr_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("NMRExperiment") %>%
    tidyr::gather("metabolite_peak", "area", -NMRExperiment) %>%
    dplyr::mutate(BaselineCorrected = TRUE),
  regions_basel_not_corr_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("NMRExperiment") %>%
    tidyr::gather("metabolite_peak", "area", -NMRExperiment) %>%
    dplyr::mutate(BaselineCorrected = FALSE)
) %>% ggplot() + geom_point(aes(x = NMRExperiment, y = area, color = metabolite_peak)) +
  facet_wrap(~BaselineCorrected)


## ------------------------------------------------------------------------
ppm_to_assign <- c(4.060960203, 3.048970634,2.405935596,0.990616851,0.986520147, 1.044258467)
identification <- nmr_identify_regions_blood (ppm_to_assign)

## ------------------------------------------------------------------------
full_spectra_matrix <- nmr_data(dataset)
print(full_spectra_matrix[1:4, 1:6]) # change it as you wish
nmr_data(dataset) <- full_spectra_matrix # Rewrite the matrix

## ------------------------------------------------------------------------
nsamp <- 12
npoints <- 20
# Create a random spectra matrix
dummy_ppm_axis <- seq(from = 0.2, to = 10, length.out = npoints)
dummy_spectra_matrix <- matrix(runif(nsamp*npoints), nrow = nsamp, ncol = npoints)
metadata <- list(external = data.frame(NMRExperiment = paste0("Sample", 1:12),
                                       DummyClass = c("a", "b"),
                                       stringsAsFactors = FALSE))
your_custom_nmr_dataset_1D <- new_nmr_dataset_1D(ppm_axis = dummy_ppm_axis,
                                                 data_1r = dummy_spectra_matrix,
                                                 metadata = metadata)
print(your_custom_nmr_dataset_1D)
plot(your_custom_nmr_dataset_1D) +
  ggtitle("Of course those random values don't make much sense...")

