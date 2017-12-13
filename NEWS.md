# NIHSnmr 1.1.0.9000

## Breaking changes

- `nmr_dataset_load` and `nmr_dataset_save` now use `readRDS` and `saveRDS` 
  instead of `load` and `save`. This is the right approach to serialize
  single R objects. If you need a script to convert previously saved datasets (created
  using `nmr_dataset_save`) please use
  `NIHSnmr:::nmr_dataset_load_old_and_save("your_old_file.RData", "your_old_file.RDS")`
  to convert the files. Sorry for the inconvenience, but the sooner we fix this
  the better.