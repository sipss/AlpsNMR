# NIHSnmr 1.1.0.9000

## Breaking changes

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

## Fixes

- Remove workaround to dplyr issue: https://github.com/tidyverse/dplyr/issues/2203
  (Sergio Oller reported and fixed the issue, dplyr-0.7.0 is fixed)

