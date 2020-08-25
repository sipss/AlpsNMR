# AlpsNMR 2.99.1 (2020-08-25)

- AlpsNMR.Rproj added to gitignore
- Modified examples to avoid create files in main package folder

# AlpsNMR 2.99.0 (2020-08-24)

- Added bootstrap and permutation method and some plots
  related to it
- Minor modifications for bioconductor submision

# AlpsNMR 2.5.9002 (2020-05-25)

- Changes to pass BiocCheck
- Added permutation test and permutation test plot to
  `nmr_data_analysis`

# AlpsNMR 2.4.9002 (2020-05-13)

- Changes to pass checks for R4

# AlpsNMR 2.3.3.9002

- NIHS_specific removed
- Tests coverage up to 30%
- Update of `save_profiling_plots`
- Add tutorial
- Remotes installation
- nmr_diagnose is deprecated. Since nmr_diagnose was only used for getting extra
  normalization information, it was been replaced with `nmr_normalize_extra_info`
  that offers a less confusing name.


# AlpsNMR 2.3.3.9001

- Add `nmr_identify_regions_cell` function
- Add documentation of `HMDB_cell`
- Vignettes updated
- New functions to apply multilevel statistics
- Update of README file

# AlpsNMR 2.3.3

- Change of `nmr_identify_regions_blood` function
- Add `nmr_identify_regions_urine` function
- Add documentation of `HMDB_urine`
- Add `computes_peak_width_ppm`function for `nmr_integrate_peak_positions`
- New `get_integration_with_metadata`
- Vignettes updated
- New functions to apply machine learning to proccessed datasets

# AlpsNMR 2.3.2

- Inclusion of baseline removal using assymetric least squares
- Change the baselineThresh to NULL so it is autodetected
- Vignettes updated including baseline removal
- Bug correction in nmr_baseline_threshold
- Elimination of package vignettes (there is an error to be solved there)
- New `nmr_identify_regions` function
- Add documentation of `HMDB_blood`
- New `files_to_rDolphin` function

# AlpsNMR 2.3.1.9000

- Rename package from NIHSnmr to AlpsNMR

# NIHSnmr 2.3.1

- Change SNR.Th value from 3 to 4 in pipeline_example.R
- Update installation instructions
- Last version form Sergio (changes not significant since 2.3.0)

# NIHSnmr 2.3.0

- Improve installation instructions under R<3.5
- nmr_peak_detection_tune_snr function added.
- Minor bug fixes

# NIHSnmr 2.2.0

- Improve installation instructions
- Clarify Add metadata vignette
- Add normalization diagnostics
- Add some data analysis helpers
- Enable parallellization for sample loading, peak detection and data analysis helpers
- Do not set negative area values to zero, to avoid biasing variances
- Add metadata from a single tidy excel function
- Add nmr_diagnose to get and set diagnostic information
- Add nmr_diagnose support to nmr_normalize
- Minor bug fixes

# NIHSnmr 2.1.0

- Documentation improvements
- nmr_dataset_peak_table object for peak detection results

# NIHSnmr 2.0.0

- Too many changes to be listed here. Check the vignette for a summary of all
the features. Use `browseVignettes("NIHSnmr")`.

# NIHSnmr 1.2.0

## Breaking changes

- Rename `injection_id` to `NMRExperiment`.

- `nmr_dataset_load` and `nmr_dataset_save` now use `readRDS` and `saveRDS` 
  instead of `load` and `save`. This is the right approach to serialize
  single R objects. If you need a script to convert previously saved datasets (created
  using `nmr_dataset_save`) please use
  `NIHSnmr:::nmr_dataset_load_old_and_save("your_old_file.RData", "your_old_file.RDS")`
  to convert the files. Sorry for the inconvenience, but the sooner we fix this
  the better.

- `filter` to select a subset of samples from an `nmr_dataset` object has
  been adapted to `dplyr >= 0.7.4`. Unless you used the `.dots` argument in
  your calls there is no need to change anything. This means we now use a tidy
  evaluation syntax for `filter`.

- `nmr_get_metadata()` returns always a data frame / tibble, even when only a single
  column is requested. It also always includes the "NMRExperiment" column.

- `nmr_dataset` object has two tables `metadata` and `metadata_ext`. The
  `metadata_ext` table includes all the metadata we add with `nmr_add_metadata` while
  `metadata` has the internal metadata (acquisition parameters, etc).
  Please use `nmr_get_metadata(nmr_dataset)` instead of `nmr_dataset$metadata`.

## Other changes

- Remove workaround to dplyr issue: https://github.com/tidyverse/dplyr/issues/2203
  (Sergio Oller reported and fixed the issue, dplyr-0.7.0 is fixed)

- The Bruker title file has quite a free format definition. A title file can
  contain lines like "Field value" or "Field value ;" or simply "value".
  The heuristics to parse the title file have been improved.
  
- Depend on tidyr 0.8.1. tidyr 0.8.0 had a bug that we reported (and for which we
  also provided a fix): https://github.com/tidyverse/tidyr/pull/419

- `nmr_get_metadata` gives a warning if the user asks for metadata columns that
  are missing.

- New `nmr_integrate_regions` function.

- `nmr_normalize` accepts `pqn` normalization.
