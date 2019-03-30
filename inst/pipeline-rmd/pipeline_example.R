# This demo script emulates a pipeline

#### Input parameters given by the user #######################################


### First node: Load Samples 

# This directory contains a sample dataset: (e.g. "/dir/subdir/dir_with_samples")
load_samples_input_dir <- system.file("dataset-demo", package = "AlpsNMR")
# A sample filtering wildcard (globbing pattern) (e.g. "*0", "*0.zip")
load_samples_glob <- "*0.zip"


### Second Node: Append metadata
# Excel file with metadata to merge:
add_metadata_excel_file <- system.file("dataset-demo", "dummy_metadata.xlsx", package = "AlpsNMR")


### Third Node: Interpolate 1D
# The range with optionally the ppm axis resolution. Typically c(0.2, 10, 2.3E-4)
interpolate1d_ppm_axis <- c(min = 0.2, max = 10, by = 2.3e-4)

### Fourth node: 
exclude_regions <-  list(water = c(4.6, 5.0), methanol = c(3.33, 3.39))

### Fifth node: First Outlier Detection
# No parameters required :-O

### Sixth node: Filter samples
samples_to_keep_conditions <- 'NMRExperiment != "40"'


#### Seventh node: Peak detection and Alignment

# Leave those as recommended for plasma.
nDivRange_ppm <- 0.1
baselineThresh <- 0
SNR.Th <- 4

maxShift_ppm <- 0.0015


#### Eighth node:  Normalization
internal_calibrant <- c(6.002, 6.015) # Can be a ppm range or NULL

#### Ninth node: Integration

peak_width_ppm <- 0.0023


#### Where the outputs of this example will be

output_dir <- tempfile(pattern = "pipeline_output_")
message("Pipeline outputs will be written at: ", output_dir)

#### Number of cores to use

num_workers <- 4

#### Pipeline begins here ######################################################
library(AlpsNMR)


#### First node: Load samples ##################################################


# Output directory:
output_dir_load <- file.path(output_dir, "01-load-samples")

plan(multiprocess, workers = num_workers)
# Load all samples ending in "0":
pipe_load_samples(samples_dir = load_samples_input_dir,
                  glob = load_samples_glob,
                  output_dir = output_dir_load)

plan(sequential)

#### Second node: Append metadata ##############################################

# Output directory:
output_dir_add_metadata <- file.path(output_dir, "02-add-metadata")

# The nmr_dataset location from the first node:
nmr_dataset_rds <- file.path(output_dir_load, "nmr_dataset.rds")

pipe_add_metadata(nmr_dataset_rds = nmr_dataset_rds,
                  output_dir = output_dir_add_metadata,
                  excel_file = add_metadata_excel_file)

#### Third node: Interpolate 1D ################################################

nmr_dataset_rds <- file.path(output_dir_add_metadata, "nmr_dataset.rds")
output_dir_interpolate1D <- file.path(output_dir, "03-interpolate-1D")
pipe_interpolate_1D(nmr_dataset_rds, interpolate1d_ppm_axis, output_dir_interpolate1D)


#### Fourth node: Exclude regions ##############################################

nmr_dataset_rds <- file.path(output_dir_interpolate1D, "nmr_dataset.rds")
output_dir_exclude <- file.path(output_dir, "04-exclude-regions")
pipe_exclude_regions(nmr_dataset_rds, exclude_regions, output_dir_exclude)


#### Fifth node: Outlier detection #############################################
nmr_dataset_rds <- file.path(output_dir_exclude, "nmr_dataset.rds")
output_dir_outlier1 <- file.path(output_dir, "05-outlier-detection")

pipe_outlier_detection(nmr_dataset_rds, output_dir_outlier1)


#### Sixth node: Filter samples ################################################

# This node is useful to filter by cohort
nmr_dataset_rds <- file.path(output_dir_outlier1, "nmr_dataset.rds")
output_dir_filter <- file.path(output_dir, "06-filter-samples")
pipe_filter_samples(nmr_dataset_rds, samples_to_keep_conditions, output_dir_filter)

#### Seventh node: Peak detection and Alignment ##################################

nmr_dataset_rds <- file.path(output_dir_filter, "nmr_dataset.rds")
output_dir_alignment <- file.path(output_dir, "07-alignment")

plan(multiprocess, workers = num_workers)
pipe_peakdet_align(nmr_dataset_rds,
                   nDivRange_ppm = nDivRange_ppm, scales = seq(1, 16, 2),
                   baselineThresh = baselineThresh, SNR.Th = SNR.Th,
                   maxShift_ppm = maxShift_ppm, acceptLostPeak = FALSE,
                   output_dir = output_dir_alignment)
plan(sequential)

#### Eighth node: PQN Normalization #######################
nmr_dataset_rds <- file.path(output_dir_alignment, "nmr_dataset.rds")
output_dir_normalization <- file.path(output_dir, "08-normalization")

pipe_normalization(nmr_dataset_rds, internal_calibrant, output_dir_normalization)


#### Ninth node: Peak integration ############################################

nmr_dataset_rds <- file.path(output_dir_normalization, "nmr_dataset.rds")
output_dir_integration <- file.path(output_dir, "09-peak-integration")

pipe_peak_integration(nmr_dataset_rds, peak_det_align_dir = output_dir_alignment,
                      peak_width_ppm = peak_width_ppm, output_dir_integration)


########################################################################
#### Finished ####
########################################################################

message("Don't forget to check out: ",  output_dir)
