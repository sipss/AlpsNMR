# Package index

## Basic functions

Functions you’ll use on a first analysis

- [`AlpsNMR`](https://sipss.github.io/AlpsNMR/reference/AlpsNMR-package.md)
  [`AlpsNMR-package`](https://sipss.github.io/AlpsNMR/reference/AlpsNMR-package.md)
  : AlpsNMR: Automated spectraL Processing System for NMR
- [`nmr_read_samples_dir()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_samples.md)
  [`nmr_read_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_samples.md)
  : Read NMR samples
- [`nmr_meta_add()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_add.md)
  [`nmr_meta_add_tidy_excel()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_add.md)
  : Add metadata to an nmr_dataset object
- [`nmr_interpolate_1D()`](https://sipss.github.io/AlpsNMR/reference/nmr_interpolate_1D.md)
  : Interpolate a set of 1D NMR Spectra
- [`nmr_exclude_region()`](https://sipss.github.io/AlpsNMR/reference/nmr_exclude_region.md)
  : Exclude region from samples
- [`plot(`*`<nmr_dataset_1D>`*`)`](https://sipss.github.io/AlpsNMR/reference/plot.nmr_dataset_1D.md)
  : Plot an nmr_dataset_1D
- [`nmr_baseline_removal()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_removal.md)
  : Baseline Removal NMR
- [`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md)
  : Peak detection for NMR
- [`nmr_align()`](https://sipss.github.io/AlpsNMR/reference/nmr_align.md)
  : Align NMR spectra
- [`nmr_normalize()`](https://sipss.github.io/AlpsNMR/reference/nmr_normalize.md)
  [`nmr_normalize_extra_info()`](https://sipss.github.io/AlpsNMR/reference/nmr_normalize.md)
  : Normalize nmr_dataset_1D samples
- [`nmr_integrate_peak_positions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_peak_positions.md)
  : Integrate peak positions
- [`nmr_data_analysis()`](https://sipss.github.io/AlpsNMR/reference/nmr_data_analysis.md)
  : Data analysis

## Metadata related functions

Functions to handle metadata

- [`nmr_meta_add()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_add.md)
  [`nmr_meta_add_tidy_excel()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_add.md)
  : Add metadata to an nmr_dataset object
- [`nmr_meta_export()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_export.md)
  : Export Metadata to an Excel file
- [`nmr_meta_get_column()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get_column.md)
  : Get a single metadata column
- [`nmr_meta_get()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_get.md)
  : Get metadata
- [`nmr_meta_groups()`](https://sipss.github.io/AlpsNMR/reference/nmr_meta_groups.md)
  : Get the names of metadata groups

## Import / Export functions

Functions to read NMR files into AlpsNMR datasets or to export to other
packages

- [`nmr_dataset_load()`](https://sipss.github.io/AlpsNMR/reference/load_and_save_functions.md)
  [`nmr_dataset_save()`](https://sipss.github.io/AlpsNMR/reference/load_and_save_functions.md)
  : nmr_dataset_load
- [`nmr_read_bruker_fid()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_bruker_fid.md)
  : Read Free Induction Decay file
- [`nmr_read_samples_dir()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_samples.md)
  [`nmr_read_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_samples.md)
  : Read NMR samples
- [`nmr_zip_bruker_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_zip_bruker_samples.md)
  : Create one zip file for each brucker sample path
- [`nmr_export_data_1r()`](https://sipss.github.io/AlpsNMR/reference/nmr_export_data_1r.md)
  : Export 1D NMR data to a CSV file
- [`files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/files_to_rDolphin.md)
  : Files to rDoplhin
- [`save_files_to_rDolphin()`](https://sipss.github.io/AlpsNMR/reference/save_files_to_rDolphin.md)
  : Save files to rDoplhin
- [`nmr_data_1r_to_SummarizedExperiment()`](https://sipss.github.io/AlpsNMR/reference/nmr_data_1r_to_SummarizedExperiment.md)
  : Export 1D NMR data to SummarizedExperiment
- [`SummarizedExperiment_to_nmr_data_1r()`](https://sipss.github.io/AlpsNMR/reference/SummarizedExperiment_to_nmr_data_1r.md)
  : Import SummarizedExperiment as 1D NMR data
- [`SummarizedExperiment_to_nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/SummarizedExperiment_to_nmr_dataset_peak_table.md)
  : Import SummarizedExperiment as mr_dataset_peak_table
- [`to_ChemoSpec()`](https://sipss.github.io/AlpsNMR/reference/to_ChemoSpec.md)
  : Convert to ChemoSpec Spectra class
- [`to_ASICS()`](https://sipss.github.io/AlpsNMR/reference/to_ASICS.md)
  : Export data for the ASICS spectral quantification library
- [`nmr_batman_options()`](https://sipss.github.io/AlpsNMR/reference/nmr_batman_options.md)
  : Batman Options helper
- [`nmr_batman_write_options()`](https://sipss.github.io/AlpsNMR/reference/nmr_batman.md)
  [`nmr_batman_export_dataset()`](https://sipss.github.io/AlpsNMR/reference/nmr_batman.md)
  [`nmr_batman_multi_data_user_hmdb()`](https://sipss.github.io/AlpsNMR/reference/nmr_batman.md)
  [`nmr_batman_multi_data_user()`](https://sipss.github.io/AlpsNMR/reference/nmr_batman.md)
  [`nmr_batman_metabolites_list()`](https://sipss.github.io/AlpsNMR/reference/nmr_batman.md)
  : Batman helpers
- [`nmr_dataset_peak_table_to_SummarizedExperiment()`](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_peak_table_to_SummarizedExperiment.md)
  : Export nmr_dataset_peak_table to SummarizedExperiment
- [`save_profiling_output()`](https://sipss.github.io/AlpsNMR/reference/save_profiling_output.md)
  : Save rDoplhin output

## Data analysis / Modelling functions

Functions for machine learning or data analysis, to be used once the
data has been preprocessed

- [`bp_kfold_VIP_analysis()`](https://sipss.github.io/AlpsNMR/reference/bp_kfold_VIP_analysis.md)
  : K-fold bootstrap and permutation over PLS-VIP
- [`bp_VIP_analysis()`](https://sipss.github.io/AlpsNMR/reference/bp_VIP_analysis.md)
  : Bootstrap and permutation over PLS-VIP
- [`models_stability_plot_bootstrap()`](https://sipss.github.io/AlpsNMR/reference/models_stability_plot_bootstrap.md)
  : Models stability plot
- [`models_stability_plot_plsda()`](https://sipss.github.io/AlpsNMR/reference/models_stability_plot_plsda.md)
  : Models stability plot
- [`new_nmr_data_analysis_method()`](https://sipss.github.io/AlpsNMR/reference/nmr_data_analysis_method.md)
  : Create method for NMR data analysis
- [`nmr_data_analysis()`](https://sipss.github.io/AlpsNMR/reference/nmr_data_analysis.md)
  : Data analysis
- [`permutation_test_model()`](https://sipss.github.io/AlpsNMR/reference/permutation_test_model.md)
  : Permutation test
- [`permutation_test_plot()`](https://sipss.github.io/AlpsNMR/reference/permutation_test_plot.md)
  : Permutation test plot
- [`plot_bootstrap_multimodel()`](https://sipss.github.io/AlpsNMR/reference/plot_bootstrap_multimodel.md)
  : Bootstrap plot predictions
- [`plot_plsda_multimodel()`](https://sipss.github.io/AlpsNMR/reference/plot_plsda_multimodel.md)
  : Multi PLDSA model plot predictions
- [`plot_plsda_samples()`](https://sipss.github.io/AlpsNMR/reference/plot_plsda_samples.md)
  : Plot PLSDA predictions
- [`plot_vip_scores()`](https://sipss.github.io/AlpsNMR/reference/plot_vip_scores.md)
  : Plot vip scores of bootstrap
- [`plsda_auroc_vip_compare()`](https://sipss.github.io/AlpsNMR/reference/plsda_auroc_vip_compare.md)
  : Compare PLSDA auroc VIP results
- [`plsda_auroc_vip_method()`](https://sipss.github.io/AlpsNMR/reference/plsda_auroc_vip_method.md)
  : Method for nmr_data_analysis (PLSDA model with AUROC and VIP
  outputs)
- [`random_subsampling()`](https://sipss.github.io/AlpsNMR/reference/random_subsampling.md)
  : Random subsampling

## nmr_dataset manipulation functions

Functions to manage and filter nmr_dataset objects

- [`filter(`*`<nmr_dataset_family>`*`)`](https://sipss.github.io/AlpsNMR/reference/filter.nmr_dataset_family.md)
  : Keep samples based on metadata column criteria
- [`format(`*`<nmr_dataset_1D>`*`)`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_1D.md)
  : format for nmr_dataset_1D
- [`format(`*`<nmr_dataset_peak_table>`*`)`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset_peak_table.md)
  : Format for nmr_dataset_peak_table
- [`format(`*`<nmr_dataset>`*`)`](https://sipss.github.io/AlpsNMR/reference/format.nmr_dataset.md)
  : Format for nmr_dataset
- [`is.nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_1D.md)
  : Object is of nmr_dataset_1D class
- [`is.nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset_peak_table.md)
  : Object is of nmr_dataset_peak_table class
- [`is.nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/is.nmr_dataset.md)
  : Object is of nmr_dataset class
- [`new_nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/new_nmr_dataset_1D.md)
  : Creates a new 1D nmr_dataset object from scratch
- [`new_nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/new_nmr_dataset_peak_table.md)
  : Creates a new nmr_dataset_peak_table object from scratch
- [`new_nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/new_nmr_dataset.md)
  : Create an nmr_dataset object
- [`nmr_data()`](https://sipss.github.io/AlpsNMR/reference/nmr_data.md)
  [`` `nmr_data<-`() ``](https://sipss.github.io/AlpsNMR/reference/nmr_data.md)
  : Set/Return the full spectra matrix
- [`nmr_dataset_1D`](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_1D.md)
  : nmr_dataset_1D (S3 class)
- [`nmr_dataset_family`](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_family.md)
  : nmr_dataset like objects (S3 classes)
- [`as.data.frame(`*`<nmr_dataset_peak_table>`*`)`](https://sipss.github.io/AlpsNMR/reference/nmr_dataset_peak_table.md)
  : nmr_dataset_peak_table (S3 class)
- [`nmr_dataset`](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md)
  : nmr_dataset (S3 class)
- [`plot_interactive()`](https://sipss.github.io/AlpsNMR/reference/plot_interactive.md)
  : Plots in WebGL
- [`plot_webgl()`](https://sipss.github.io/AlpsNMR/reference/plot_webgl.md)
  : Plot a dataset into a HTML file
- [`plot(`*`<nmr_dataset_1D>`*`)`](https://sipss.github.io/AlpsNMR/reference/plot.nmr_dataset_1D.md)
  : Plot an nmr_dataset_1D
- [`print(`*`<nmr_dataset_1D>`*`)`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_1D.md)
  : print for nmr_dataset_1D
- [`tidy(`*`<nmr_dataset_1D>`*`)`](https://sipss.github.io/AlpsNMR/reference/tidy.nmr_dataset_1D.md)
  : Get a tidy data frame from nmr_data object
- [`print(`*`<nmr_dataset_peak_table>`*`)`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset_peak_table.md)
  : print for nmr_dataset_peak_table
- [`print(`*`<nmr_dataset>`*`)`](https://sipss.github.io/AlpsNMR/reference/print.nmr_dataset.md)
  : Print for nmr_dataset
- [`` `[`( ``*`<nmr_dataset_1D>`*`)`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset_1D.md)
  : Extract parts of an nmr_dataset_1D
- [`` `[`( ``*`<nmr_dataset_peak_table>`*`)`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset_peak_table.md)
  : Extract parts of an nmr_dataset_peak_table
- [`` `[`( ``*`<nmr_dataset>`*`)`](https://sipss.github.io/AlpsNMR/reference/sub-.nmr_dataset.md)
  : Extract parts of an nmr_dataset
- [`validate_nmr_dataset_family()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset_family.md)
  : Validate nmr_dataset_family objects
- [`validate_nmr_dataset_peak_table()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset_peak_table.md)
  : Validate nmr_dataset_peak_table objects
- [`validate_nmr_dataset()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset.md)
  [`validate_nmr_dataset_1D()`](https://sipss.github.io/AlpsNMR/reference/validate_nmr_dataset.md)
  : Validate nmr_dataset objects

## Metabolite identification

Functions and tables related to metabolite identification

- [`HMDB_blood`](https://sipss.github.io/AlpsNMR/reference/HMDB_blood.md)
  : The Human Metabolome DataBase multiplet table: blood metabolites
  normally found in NMR-based metabolomics
- [`HMDB_cell`](https://sipss.github.io/AlpsNMR/reference/HMDB_cell.md)
  : The Human Metabolome DataBase multiplet table: cell metabolites
  normally found in NMR-based metabolomics
- [`HMDB_urine`](https://sipss.github.io/AlpsNMR/reference/HMDB_urine.md)
  : The Human Metabolome DataBase multiplet table: urine metabolites
  normally found in NMR-based metabolomics
- [`hmdb`](https://sipss.github.io/AlpsNMR/reference/hmdb.md) : The
  Human Metabolome DataBase multiplet table
- [`Parameters_blood`](https://sipss.github.io/AlpsNMR/reference/Parameters_blood.md)
  : to rDolphin
- [`Parameters_cell`](https://sipss.github.io/AlpsNMR/reference/Parameters_cell.md)
  : Parameters for cell samples profiling
- [`Parameters_urine`](https://sipss.github.io/AlpsNMR/reference/Parameters_urine.md)
  : Parameters for urine samples profiling
- [`ROI_blood`](https://sipss.github.io/AlpsNMR/reference/ROI_blood.md)
  : ROIs for blood (plasma/serum) samples
- [`ROI_cell`](https://sipss.github.io/AlpsNMR/reference/ROI_cell.md) :
  ROIs for cell samples
- [`ROI_urine`](https://sipss.github.io/AlpsNMR/reference/ROI_urine.md)
  : ROIs for urine samples
- [`nmr_identify_regions_blood()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_blood.md)
  : NMR peak identification (plasma/serum samples)
- [`nmr_identify_regions_cell()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_cell.md)
  : NMR peak identification (cell samples)
- [`nmr_identify_regions_urine()`](https://sipss.github.io/AlpsNMR/reference/nmr_identify_regions_urine.md)
  : NMR peak identification (urine samples)

## Alignment functions

Functions to align spectra

- [`nmr_align_find_ref()`](https://sipss.github.io/AlpsNMR/reference/nmr_align_find_ref.md)
  : Find alignment reference
- [`nmr_align()`](https://sipss.github.io/AlpsNMR/reference/nmr_align.md)
  : Align NMR spectra

## Baseline correction

Functions to remove or estimate the baseline

- [`nmr_baseline_removal()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_removal.md)
  : Baseline Removal NMR

- [`nmr_baseline_estimation()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_estimation.md)
  :

  Estimate the baseline on an nmr_dataset_1D object, using
  [baseline::baseline.als](https://rdrr.io/pkg/baseline/man/baseline.als.html).

## Peak detection

Functions to detect peaks and visualize the resulting peak lists

- [`nmr_detect_peaks_plot_overview()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot_overview.md)
  : Overview of the peak detection results

- [`nmr_detect_peaks_plot_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot_peaks.md)
  : Plot multiple peaks from a peak list

- [`nmr_detect_peaks_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_plot.md)
  : Plot peak detection results

- [`nmr_detect_peaks_tune_snr()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks_tune_snr.md)
  : Diagnose SNR threshold in peak detection

- [`nmr_detect_peaks()`](https://sipss.github.io/AlpsNMR/reference/nmr_detect_peaks.md)
  : Peak detection for NMR

- [`nmr_baseline_threshold_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold_plot.md)
  : Plot the baseline thresholds

- [`nmr_baseline_threshold()`](https://sipss.github.io/AlpsNMR/reference/nmr_baseline_threshold.md)
  : Threshold estimation for peak detection

- [`Peak_detection`](https://sipss.github.io/AlpsNMR/reference/Peak_detection.md)
  : Peak detection for NMR

- [`peaklist_accept_peaks()`](https://sipss.github.io/AlpsNMR/reference/peaklist_accept_peaks.md)
  :

  Peak list: Create an `accepted` column based on some criteria

- [`peaklist_fit_lorentzians()`](https://sipss.github.io/AlpsNMR/reference/peaklist_fit_lorentzians.md)
  : Fit lorentzians to each peak to estimate areas

## Peak clustering / matching / Peak tables

Functions to build peak tables from peak lists

- [`nmr_peak_clustering()`](https://sipss.github.io/AlpsNMR/reference/nmr_peak_clustering.md)
  : Peak clustering
- [`nmr_peak_clustering_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_peak_clustering_plot.md)
  : Plot clustering results
- [`nmr_build_peak_table()`](https://sipss.github.io/AlpsNMR/reference/nmr_build_peak_table.md)
  : Build a peak table from the clustered peak list
- [`nmr_get_peak_distances()`](https://sipss.github.io/AlpsNMR/reference/nmr_get_peak_distances.md)
  : Compute peak to peak distances
- [`nmr_integrate_regions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_regions.md)
  : Integrate regions
- [`nmr_integrate_peak_positions()`](https://sipss.github.io/AlpsNMR/reference/nmr_integrate_peak_positions.md)
  : Integrate peak positions

## Outlier detection

Functions to detect and report outliers in the dataset

- [`nmr_pca_build_model()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_build_model.md)
  : Build a PCA on for an nmr_dataset
- [`nmr_pca_outliers_filter()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_filter.md)
  : Exclude outliers
- [`nmr_pca_outliers_plot()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_plot.md)
  : Plot for outlier detection diagnostic
- [`nmr_pca_outliers_robust()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers_robust.md)
  : Outlier detection through robust PCA
- [`nmr_pca_outliers()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_outliers.md)
  : Compute PCA residuals and score distance for each sample
- [`nmr_pca_plot_variance()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_plots.md)
  [`nmr_pca_scoreplot()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_plots.md)
  [`nmr_pca_loadingplot()`](https://sipss.github.io/AlpsNMR/reference/nmr_pca_plots.md)
  : Plotting functions for PCA

## Other functions

Other, not catalogued, AlpsNMR functions

- [`download_MTBLS242()`](https://sipss.github.io/AlpsNMR/reference/download_MTBLS242.md)
  : Download MTBLS242

- [`file_lister()`](https://sipss.github.io/AlpsNMR/reference/file_lister.md)
  : NMR file lister

- [`get_integration_with_metadata()`](https://sipss.github.io/AlpsNMR/reference/get_integration_with_metadata.md)
  :

  Get integrals with metadata from `integrate peak positions`

- [`pipe_load_samples()`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md)
  [`pipe_add_metadata()`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md)
  [`pipe_interpolate_1D()`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md)
  [`pipe_exclude_regions()`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md)
  [`pipe_outlier_detection()`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md)
  [`pipe_filter_samples()`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md)
  [`pipe_peakdet_align()`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md)
  [`pipe_peak_integration()`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md)
  [`pipe_normalization()`](https://sipss.github.io/AlpsNMR/reference/Pipelines.md)
  : Pipelines

- [`nmr_ppm_resolution()`](https://sipss.github.io/AlpsNMR/reference/nmr_ppm_resolution.md)
  : PPM resolution of the spectra

- [`ppm_resolution()`](https://sipss.github.io/AlpsNMR/reference/ppm_resolution.md)
  : Unlisted PPM resolution

- [`nmr_autophase()`](https://sipss.github.io/AlpsNMR/reference/nmr_autophase.md)
  : Rephase 1D NMR data
