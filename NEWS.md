# AlpsNMR 3.5.1 (2022-04-07)

- `plot_interactive` now accepts an `overwrite` argument to avoid asking the user
  interactively
- Improve `nmr_detect_peaks_tune_snr` to tune the SNR threshold with the right
  other parameters
- Documentation improvements
  * Split Peak_detection page into smaller and more specific pages
- Let the user choose how code is parallellized, as suggested by BiocParallel
  documentation.
- Replace furr/future parallellization loops with BiocParallel.
  Provides a warning in case a future::plan() has been set.
- Demote Imports to Suggests: SummarizedExpriment, S4Vectors, ggrepel, GGally
- Remove dependencies: tidyselect, assertthat, plyr, furrr
- Add `download_MTBLS242()` function to help download the data for the tutorial
- Skip mixOmics test if affected by https://github.com/mixOmicsTeam/mixOmics/pull/199
- Fix auto setting of the baseline threshold for the peak detection


# AlpsNMR 3.3.4 (2021-09-16)

- Fix issue with PCA plots not working as expected
- Ensure NMRExperiment names are not duplicated in a dataset (closes #44)
- Fix issue with some title file formatting in Bruker samples (closes #46)
- Export groups in to_ChemoSpec
- License since AlpsNMR was released has alwayd been MIT
  as stated in the bioinformatics paper

# AlpsNMR 3.1.5 (2021-3-31)
- Removed warning about future_options deprecation

# AlpsNMR 3.1.4 (2021-1-20)
- bug fix loading bruker files

# AlpsNMR 3.1.3 (2020-11-19)

- Added instructions to follow a longer tutorial
- nmr_pca_outliers_plot modified to show names in all
  boundaries of the plot

# AlpsNMR 3.1.2 (2020-11-04)

- Bug fix related with Bioconductor Renviron variable
_R_CHECK_LENGTH_1_CONDITION_

# AlpsNMR 3.1.1 (2020-10-30)

- Modified order of author list

# AlpsNMR 3.1.0 (2020-10-22)

- Package accepted in bioconductor

# AlpsNMR 2.99.93 (2020-10-22)

- Héctor removed as maintainer to ensure a single point of contact

# AlpsNMR 2.99.92 (2020-10-22)

- Héctor added as maintainer

# AlpsNMR 2.99.91 (2020-10-22)

- test changed

# AlpsNMR 2.99.9 (2020-10-22)

- Added Héctor ass author
- Bug fix in nmr_read_bruker_fid

# AlpsNMR 2.99.8 (2020-10-22)

- Deletion of tutorial and download_MTBLS242_demo.R

# AlpsNMR 2.99.7 (2020-10-19)

- Bugs in import/export functions to SummarizedExperiment solved

# AlpsNMR 2.99.6 (2020-10-19)

- Added import/export options form nmr_dataset_1D to SummarizedExperiment
- Added import/export options form nmr_dataset_peak_table to SummarizedExperiment

# AlpsNMR 2.99.5 (2020-10-14)

- Bug in bp_kfold_VIP_analysis solved
- Several packages moved from import to depends
- Reexport of some functions removed
- to_rDolphin_blood code reorganized
- Typos removed from tutorial
- norm_pqn_diagnostic$norm_factor used in tutorial instead of plot it
- Parallel changed for BiocParallel

# AlpsNMR 2.99.4 (2020-09-28)

- Warning in plot_interactive function added
- Suppressed other warnings of plot_interactive function

# AlpsNMR 2.99.3 (2020-09-21)

- sapply calls changed for vapply
- Bioconductor installation instructions included
- MIT license removed
- LazyData: TRUE removed
- Excessive print statements removed from vignettes
- sessionInfo() added to end of vignettes
- Created inst/script directoy to describe inst/extdata source and
  creation #TODO falta rellenar el archivo
- Commented out code removed

# AlpsNMR 2.99.2 (2020-08-26)

- AlpsNMR.Rproj removed from git repository
- Reduced demo dataset to avoid package size > 5 MB
- Modified introduction to alpsnmr vignette and some tests
  to work with reduced demo dataset

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
