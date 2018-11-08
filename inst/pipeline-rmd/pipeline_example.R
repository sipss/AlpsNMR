# This demo script emulates a pipeline

library(NIHSnmr)
# This directory contains a sample dataset:
directory_with_samples <- system.file("dataset-demo", package = "NIHSnmr")
message("Sample dataset found at: ", directory_with_samples)
output_dir <- tempdir()
message("Pipeline outputs will be written at: ", output_dir)

############################################################################
## First node: Load samples
############################################################################

# Output directory:
output_dir_load <- file.path(output_dir, "01-load-samples")

# Load all samples ending in "0":
pipe_load_samples(samples_dir = directory_with_samples,
                  output_dir = output_dir_load,
                  glob = "*0.zip")

############################################################################
## Second node: Append metadata
############################################################################
# Here the user provides an Excel file as described at ?pipe_add_metadata

# Output directory:
output_dir_add_metadata <- file.path(output_dir, "02-add-metadata")

# Excel file with metadata to merge:
excel_file <- system.file("dataset-demo", "dummy_metadata.xlsx", package = "NIHSnmr")

# The nmr_dataset location from the first node:
nmr_dataset_rds <- file.path(output_dir_load, "nmr_dataset.rds")

pipe_add_metadata(nmr_dataset_rds = nmr_dataset_rds,
                  output_dir = output_dir_add_metadata,
                  excel_file = excel_file)

#########################################################################
## Third node: Interpolate 1D
#########################################################################
nmr_dataset_rds <- file.path(output_dir_add_metadata, "nmr_dataset.rds")
ppm_axis <- c(min = 0.2, max = 10, by = 8e-4)
output_dir_interpolate1D <- file.path(output_dir, "03-interpolate-1D")
pipe_interpolate_1D(nmr_dataset_rds, ppm_axis, output_dir_interpolate1D)

#########################################################################
## Fourth node: Exclude regions
#########################################################################
nmr_dataset_rds <- file.path(output_dir_interpolate1D, "nmr_dataset.rds")
exclude_regions <-  list(water = c(4.6, 5.0), methanol = c(3.33, 3.39))
output_dir_exclude <- file.path(output_dir, "04-exclude-regions")
pipe_exclude_regions(nmr_dataset_rds, exclude_regions, output_dir_exclude)



message("Don't forget to check out: ",  output_dir)

