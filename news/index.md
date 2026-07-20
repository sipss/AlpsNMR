# Changelog

## AlpsNMR (development version)

### Breaking changes

- [`plot_vip_scores()`](https://sipss.github.io/AlpsNMR/reference/plot_vip_scores.md):
  `nbootstrap` argument renamed to `n_samples`, and its threshold line
  now uses the correct degrees of freedom (`df = n_samples - 1`),
  matching the selection
  [`bp_VIP_analysis()`](https://sipss.github.io/AlpsNMR/reference/bp_VIP_analysis.md)
  actually performs.
- [`nmr_read_bruker_fid()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_bruker_fid.md):
  `endian` argument removed (byte order is now auto-detected from
  `acqus`), and it now returns a data frame with `time_s`/`fid_complex`
  columns instead of a raw numeric vector. Also fixes silent truncation
  of the FID to half its length.
- [`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md)
  /
  [`nmr_baseline_threshold_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold_plot.md):
  `range_without_peaks` / `chemshift_range` no longer default to
  `c(9.5, 10)` ppm; must now be given explicitly.
- [`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md):
  `range_without_peaks` no longer defaults to `c(9.5, 10)` ppm either;
  either it or `baselineThresh` must now be given, or the call aborts
  with a clear message instead of failing deep inside with a confusing
  one ([\#66](https://github.com/sipss/AlpsNMR/issues/66)).
- [`download_MTBLS242()`](https://sipss.github.io/AlpsNMR/reference/download_MTBLS242.md):
  `keep_only_preop_and_3months` (logical) replaced by `timepoints` (a
  character vector of `TimePoint` values to keep, or `NULL` for every
  timepoint), so any two (or more) of the study’s five timepoints can be
  selected, not just preop and 12 months.
  `timepoints = c("preop", "12 months after surgery")` is the new
  default and reproduces the old `keep_only_preop_and_3months = TRUE`
  behavior.
- [`nmr_baseline_threshold_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold_plot.md):
  `NMRExperiment = NULL` now means “every sample” (paginated via new
  `nrow`/`ncol`/`page` arguments) instead of silently subsampling to 10
  random samples when there were more than 20. `"all"` remains a synonym
  for `NULL`. `nrow`/`ncol` default to a snug grid for small sample
  counts, or a fixed 3x3 (paginate with `page`) for 7 or more.

### New features

- [`nmr_baseline_threshold_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold_plot.md):
  new `nrow`/`ncol`/`page` arguments paginate the per-sample facets
  instead of cramming every sample onto one illegible page.

### Bug fixes

- [`bp_kfold_VIP_analysis()`](https://sipss.github.io/AlpsNMR/reference/bp_kfold_VIP_analysis.md):
  fixed fold partitioning, which assigned samples to folds with a
  deterministic `x %% k` split instead of the intended random shuffle.
- [`bp_VIP_analysis()`](https://sipss.github.io/AlpsNMR/reference/bp_VIP_analysis.md):
  a degenerate (single-class) bootstrap resample is now redrawn instead
  of being patched with a fixed, non-random replacement, which biased a
  fraction of the bootstrap replicates.
- `read_bruker_param()`: fixed a wrong regex capture-group index causing
  `subscript out of bounds` on some Bruker parameter files.
- `parse_title_file()`: fixed a
  [`gsub()`](https://rdrr.io/r/base/grep.html) call missing
  `perl = TRUE`, which left trailing whitespace untrimmed from Bruker
  pdata title fields.
- `choose_best_nlv()`: fixed a key-name mismatch that made
  `diagnostic_plot`, `diagnostic_box_plot`, and `model_performances`
  always `NULL` for models built with
  [`plsda_auroc_vip_method()`](https://sipss.github.io/AlpsNMR/reference/plsda_auroc_vip_method.md).
- [`models_stability_plot_plsda()`](https://sipss.github.io/AlpsNMR/reference/models_stability_plot_plsda.md)
  /
  [`models_stability_plot_bootstrap()`](https://sipss.github.io/AlpsNMR/reference/models_stability_plot_bootstrap.md)
  /
  [`plot_bootstrap_multimodel()`](https://sipss.github.io/AlpsNMR/reference/plot_bootstrap_multimodel.md):
  fixed two ggplot2 arguments deprecated since ggplot2 3.3.4/3.4.0,
  which warned on every call.
- `create_sample_names()`: disambiguates samples sharing a leaf
  directory name (e.g. Bruker EXPNO `10`) by stripping the full common
  path prefix, instead of only one parent level. Collisions more than
  one level deep now get readable names instead of `vctrs`-generated
  `...N` suffixes ([\#62](https://github.com/sipss/AlpsNMR/issues/62)).
- [`tidy.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/tidy.nmr_dataset_1D.md)
  (and thus [`plot()`](https://rdrr.io/r/graphics/plot.default.html) /
  [`nmr_baseline_threshold_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold_plot.md)):
  an unknown `NMRExperiment` value now warns and is excluded (if some
  values are valid) or errors (if none are), instead of silently
  returning rows with `NA` intensities
  ([\#69](https://github.com/sipss/AlpsNMR/issues/69)).

### Other changes

- [`download_MTBLS242()`](https://sipss.github.io/AlpsNMR/reference/download_MTBLS242.md):
  validate downloaded files against MetaboLights’ published SHA-256
  checksums for the MTBLS242 dataset
  ([\#72](https://github.com/sipss/AlpsNMR/issues/72)).
- Bumped several dependency version floors to roughly their versions
  from a year ago. Packages with a recent major release are pinned to
  the last minor of the previous major instead, to avoid forcing an
  upgrade: `ggplot2 (>= 3.5.2)`, `fs (>= 1.6.7)`, `curl (>= 6.4.0)`,
  `zip (>= 2.3.3)`, `progressr (>= 0.19.0)`.

## AlpsNMR 4.11.1 (2025-09-24)

- Compatibility with ggplot2-4.0.

## AlpsNMR 4.7.2 (2024-08-10)

- Disable nested parallellization in
  [`nmr_detect_peaks_tune_snr()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_tune_snr.md).

## AlpsNMR 4.7.1 (2024-06-02)

- Added
  [`nmr_autophase()`](https://sipss.github.io/AlpsNMR/reference/nmr_autophase.md)
  for automated phase correction using the NMRphasing package
  ([\#68](https://github.com/sipss/AlpsNMR/issues/68)).
- Added `to_ASICS` function to export dataset for ASICS quantification
  ([\#68](https://github.com/sipss/AlpsNMR/issues/68)).

## AlpsNMR 4.1.6 (2023-02-16)

- Download improvements:
  - Progress bar
  - Detect user interruptions
  - Sleep 3 seconds between retrying failed downloads

## AlpsNMR 4.1.5 (2023-02-10)

- Replace deprecated dplyr::select() calls.
- Remove workaround for mixOmics bug, bump mixOmics dependency
- Simplify implementation (same algorithm) for determining the optimal
  number of latent variables in the plsda models.
- Bump dplyr dependency version.

## AlpsNMR 4.1.4 (2022-11-08)

- Disable nested parallelization
- Update workaround Biocparallel bpmapply
- Revert bpstop to sleep

## AlpsNMR 4.1.3 (2022-11-07)

- Closer to the fix

## AlpsNMR 4.1.2 (2022-11-04)

- Try a more robust fix on palomino4 (bpstop() instead sleep)
- Use register() in an example to avoid further breakage on palomino
- Workaround performance issues on BiocParallel::bpmapply()
  (<https://github.com/Bioconductor/BiocParallel/pull/228>)

## AlpsNMR 4.1.1 (2022-11-02)

- Remove archive dependency
- Try fixing build on palomino4, due to race condition in R CMD check

## AlpsNMR 3.99.7 (2022-10-27)

### Minor changes

- Fix build issue on palomino4, simplifying helper function

## AlpsNMR 3.99.6 (2022-10-26)

### Minor changes

- When saving, normalize extra information is saved as well.
- Updated downsampled demo data for examples.
- More robust nmr_baseline_threshold()
- Faster examples

## AlpsNMR 3.99.5 (2022-10-26)

### Minor changes

- Remove call to deprecated
  [`ggplot2::qplot()`](https://ggplot2.tidyverse.org/reference/qplot.html)

## AlpsNMR 3.99.4 (2022-10-19)

### Major changes

- Improved the
  [`download_MTBLS242()`](https://sipss.github.io/AlpsNMR/reference/download_MTBLS242.md)
  function, allowing to either download the parts of MTBLS242 needed for
  the tutorial or the whole dataset, which may be nice to have if you
  want to play beyond the tutorial.

- When reading a Bruker sample from a zip file, you now can specify in
  the file name the zip subdirectory. For instance,
  “/path/to/sample.zip!/sample/3”, when `sample.zip` contains a folder
  named `sample` with a subfolder named `3` that includes the sample
  data you want to actually read.

### Minor changes

- Remove Bioconductor Build System workaround, since
  <https://github.com/Bioconductor/BBS/issues/220> was fixed.

## AlpsNMR 3.99.3 (2022-10-17)

- Add libarchive as a SystemRequirement to workaround a limitation of
  the Bioconductor build system (BBS), that can’t pick system
  requirements recursively. Thanks to Jennifer Wokaty for checking the
  BBS and providing this suggestion.

## AlpsNMR 3.99.2 (2022-10-14)

### Breaking changes

- Set `fix_baseline = FALSE` in
  [`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)
  as default. The former `TRUE` approach here did not make much sense if
  peak boundaries were not perfectly established.

### Major changes

- Baseline estimation: We now offer
  [`nmr_baseline_estimation()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_estimation.md)
  besides
  [`nmr_baseline_removal()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_removal.md).
  The estimation function computes the baseline and saves it instead of
  subtracting it from the signal. This is a better approach because it
  lets each step of the pipeline decide whether it makes sense to
  subtract the baseline or not. The
  [`nmr_baseline_removal()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_removal.md)
  is for now still available, but it will be deprecated in a future
  version.

- For the `baselineThresh` argument in
  [`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md)
  we now suggest using
  `nmr_baseline_threshold(dataset, method = "median3mad")`. This is more
  robust than the former (but still the default) method.

- Peak detection and integration: We want to approach the peak
  detection, clustering an integration in a different way. While the old
  pipeline still works as expected, we have introduced new arguments to
  peak detection, with backwards compatible defaults and a peak
  clustering function. We still provide the vignette with the former
  workflow, because it is still relevant but we may deprecate it in a
  future version, once we are confident the changes we are making are
  robust across several datasets.

- Parallellization: We are switching from the `future` package to
  `BiocParallel`, to better integrate in the Bioconductor ecosystem. In
  this version, if you use a different future plan you may get a warning
  to switch to BiocParallel. In a future version we will remove our
  dependency with the (awesome) `future` package.

### Minor changes

- You can now set experiment names (NMRExperiment) with
  `names(dataset) <- c("Sample1", "Sample2")`.
- You can now pass a named vector with the sample names to the
  nmr_read_samples function. The names will be used as the sample names.
- Peak detection has a more robust baseline threshold estimation
- Peak detection estimates the baseline threshold on each sample
  individually. The threshold is calculated using only the sample where
  we are currently detecting the peaks.
- Peak detection includes a simple but effective lorentzian fitting (for
  area and width estimation)
- Add functions to evaluate the quality of the peak detection using
  plots
- More fine grained interpolation axis if `axis = NULL` is given in
  [`nmr_interpolate_1D()`](https://sipss.github.io/AlpsNMR/reference/nmr_interpolate_1D.md)
- Save list of excluded regions in the `nmr_dataset` object.
- Drop MassSpecWavelet workaround on partial argument matching since it
  was fixed upstream
- Documentation: Start providing verbose messages with tips in functions
- Remove unused deprecated imports from the `future` package
  ([\#65](https://github.com/sipss/AlpsNMR/issues/65), thanks to
  [@HenrikBengtsson](https://github.com/HenrikBengtsson))
- Add URL and BugReports to the DESCRIPTION
  ([\#64](https://github.com/sipss/AlpsNMR/issues/64), thanks to
  [@HenrikBengtsson](https://github.com/HenrikBengtsson))
- Reading bruker samples is now a bit more robust and gives detailed
  tracebacks in case of error.

## AlpsNMR 3.5.1 (2022-04-07)

- `plot_interactive` now accepts an `overwrite` argument to avoid asking
  the user interactively
- Improve `nmr_detect_peaks_tune_snr` to tune the SNR threshold with the
  right other parameters
- Documentation improvements
  - Split Peak_detection page into smaller and more specific pages
- Let the user choose how code is parallellized, as suggested by
  BiocParallel documentation.
- Replace furr/future parallellization loops with BiocParallel. Provides
  a warning in case a future::plan() has been set.
- Demote Imports to Suggests: SummarizedExpriment, S4Vectors, ggrepel,
  GGally
- Remove dependencies: tidyselect, assertthat, plyr, furrr
- Add
  [`download_MTBLS242()`](https://sipss.github.io/AlpsNMR/reference/download_MTBLS242.md)
  function to help download the data for the tutorial
- Skip mixOmics test if affected by
  <https://github.com/mixOmicsTeam/mixOmics/pull/199>
- Fix auto setting of the baseline threshold for the peak detection

## AlpsNMR 3.3.4 (2021-09-16)

- Fix issue with PCA plots not working as expected
- Ensure NMRExperiment names are not duplicated in a dataset (closes
  [\#44](https://github.com/sipss/AlpsNMR/issues/44))
- Fix issue with some title file formatting in Bruker samples (closes
  [\#46](https://github.com/sipss/AlpsNMR/issues/46))
- Export groups in to_ChemoSpec
- License since AlpsNMR was released has alwayd been MIT as stated in
  the bioinformatics paper

## AlpsNMR 3.1.5 (2021-3-31)

- Removed warning about future_options deprecation

## AlpsNMR 3.1.4 (2021-1-20)

- bug fix loading bruker files

## AlpsNMR 3.1.3 (2020-11-19)

- Added instructions to follow a longer tutorial
- nmr_pca_outliers_plot modified to show names in all boundaries of the
  plot

## AlpsNMR 3.1.2 (2020-11-04)

- Bug fix related with Bioconductor Renviron variable
  *R_CHECK_LENGTH_1_CONDITION*

## AlpsNMR 3.1.1 (2020-10-30)

- Modified order of author list

## AlpsNMR 3.1.0 (2020-10-22)

- Package accepted in bioconductor

## AlpsNMR 2.99.93 (2020-10-22)

- Héctor removed as maintainer to ensure a single point of contact

## AlpsNMR 2.99.92 (2020-10-22)

- Héctor added as maintainer

## AlpsNMR 2.99.91 (2020-10-22)

- test changed

## AlpsNMR 2.99.9 (2020-10-22)

- Added Héctor ass author
- Bug fix in nmr_read_bruker_fid

## AlpsNMR 2.99.8 (2020-10-22)

- Deletion of tutorial and download_MTBLS242_demo.R

## AlpsNMR 2.99.7 (2020-10-19)

- Bugs in import/export functions to SummarizedExperiment solved

## AlpsNMR 2.99.6 (2020-10-19)

- Added import/export options form nmr_dataset_1D to
  SummarizedExperiment
- Added import/export options form nmr_dataset_peak_table to
  SummarizedExperiment

## AlpsNMR 2.99.5 (2020-10-14)

- Bug in bp_kfold_VIP_analysis solved
- Several packages moved from import to depends
- Reexport of some functions removed
- to_rDolphin_blood code reorganized
- Typos removed from tutorial
- norm_pqn_diagnostic\$norm_factor used in tutorial instead of plot it
- Parallel changed for BiocParallel

## AlpsNMR 2.99.4 (2020-09-28)

- Warning in plot_interactive function added
- Suppressed other warnings of plot_interactive function

## AlpsNMR 2.99.3 (2020-09-21)

- sapply calls changed for vapply
- Bioconductor installation instructions included
- MIT license removed
- LazyData: TRUE removed
- Excessive print statements removed from vignettes
- sessionInfo() added to end of vignettes
- Created inst/script directoy to describe inst/extdata source and
  creation \#TODO falta rellenar el archivo
- Commented out code removed

## AlpsNMR 2.99.2 (2020-08-26)

- AlpsNMR.Rproj removed from git repository
- Reduced demo dataset to avoid package size \> 5 MB
- Modified introduction to alpsnmr vignette and some tests to work with
  reduced demo dataset

## AlpsNMR 2.99.1 (2020-08-25)

- AlpsNMR.Rproj added to gitignore
- Modified examples to avoid create files in main package folder

## AlpsNMR 2.99.0 (2020-08-24)

- Added bootstrap and permutation method and some plots related to it
- Minor modifications for bioconductor submision

## AlpsNMR 2.5.9002 (2020-05-25)

- Changes to pass BiocCheck
- Added permutation test and permutation test plot to
  `nmr_data_analysis`

## AlpsNMR 2.4.9002 (2020-05-13)

- Changes to pass checks for R4

## AlpsNMR 2.3.3.9002

- NIHS_specific removed
- Tests coverage up to 30%
- Update of `save_profiling_plots`
- Add tutorial
- Remotes installation
- nmr_diagnose is deprecated. Since nmr_diagnose was only used for
  getting extra normalization information, it was been replaced with
  `nmr_normalize_extra_info` that offers a less confusing name.

## AlpsNMR 2.3.3.9001

- Add `nmr_identify_regions_cell` function
- Add documentation of `HMDB_cell`
- Vignettes updated
- New functions to apply multilevel statistics
- Update of README file

## AlpsNMR 2.3.3

- Change of `nmr_identify_regions_blood` function
- Add `nmr_identify_regions_urine` function
- Add documentation of `HMDB_urine`
- Add `computes_peak_width_ppm`function for
  `nmr_integrate_peak_positions`
- New `get_integration_with_metadata`
- Vignettes updated
- New functions to apply machine learning to proccessed datasets

## AlpsNMR 2.3.2

- Inclusion of baseline removal using assymetric least squares
- Change the baselineThresh to NULL so it is autodetected
- Vignettes updated including baseline removal
- Bug correction in nmr_baseline_threshold
- Elimination of package vignettes (there is an error to be solved
  there)
- New `nmr_identify_regions` function
- Add documentation of `HMDB_blood`
- New `files_to_rDolphin` function

## AlpsNMR 2.3.1.9000

- Rename package from NIHSnmr to AlpsNMR
