# This demo script emulates a pipeline

library(NIHSnmr)
# This directory contains a sample dataset:
directory_with_samples <- system.file("dataset-demo", package = "NIHSnmr")
message("Sample dataset found at: ", directory_with_samples)
output_dir <- tempdir()
message("Pipeline outputs will be written at: ", output_dir)


## First node: Load samples

# Output directory:
output_dir_load <- file.path(output_dir, "01-load-samples")

# Load all samples ending in "0":
pipe_load_samples(samples_dir = directory_with_samples,
                  output_dir = output_dir_load,
                  glob = "*0.zip")

## Second node: Append metadata
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
