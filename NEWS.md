# NIHSnmr 1.2.0.9000

- Allow reversed ranges in `nmr_exclude_region`.
- Fix bug in `nmr_pca_plot_variance` when `ncomp` was fixed when building the model
- Add `nmr_zip_bruker_samples` to zip Bruker samples given their paths
- `nmr_integrate_regions` removes a baseline by default. Use `fix_baseline = FALSE` to use previous behaviour.

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
