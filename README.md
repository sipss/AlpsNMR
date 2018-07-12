NIHSnmr
=======

`NIHSnmr` is an R package that can:

-   Load NMR samples from Bruker instrumentation (both 1-D and 2-D,
    partial implementation for higher dimensions) (`nmr_read_samples`)
-   Plot 1-D NMR datasets or 2-D samples. (`plot`)
-   Perform very basic preprocessing tasks:
-   Interpolate NMR samples so they share a common axis
    (`nmr_interpolate`)
-   Exclude regions of the NMR spectra that correspond to solvents or
    undesired analytes (`nmr_exclude_region`)
-   Normalize the samples (if given a specific list of areas)
    (`nmr_normalize`)

Additionally the package feature features:

-   Support for reading only the sample metadata
    (`nmr_read_samples(..., metadata_only = TRUE)`)
-   Support for reading all samples in a directory (optionally of a
    specific pulse sequence) (`nmr_read_samples_dir`)
-   Support for iRods:
-   Search for NMR samples (`nmr_irods_search`)
-   Download and read samples (`nmr_read_samples_irods`)
-   Extract the metadata that irods requires from an `nmr_dataset`
    (`nmr_get_irods_meta`)
-   Prepare zip files for uploading a dataset to irods
    (`nmr_prepare_zip_files`)
-   Push a dataset to irods (`nmr_push_to_irods`)
-   Read exported excel files from sLims (`read_exported_slims`)
-   Load samples and view them interactively (`NIHSnmr_interactive`)

Installation
------------

Depending on your access permissions, one of this options will be
suitable:

### At IBEC

1. Go to: https://gitnose.ibec.local/compmetabolomics/NIHSnmr

2. Download zip file: https://gitnose.ibec.local/compmetabolomics/NIHSnmr/-/archive/master/NIHSnmr-master.zip

3. `devtools::install_local("C:/blahblahbla/NIHSnmr-master.zip")`

### At NIHS

    # Zip file:
    devtools::install_local("C:/blahblahbla/NIHSnmr.zip")

    # Rstudio server instance
    devtools::install_git('/nihs/Bioinformatics_home/rdollermse/git/NIHSnmr/')

    # Visual studio git repository (https):
    devtools::install_git('https://nestle-eur.visualstudio.com/NIHS_NEMO2/_git/NIHSnmr')

    # Visual studio git repository (ssh):
    devtools::install_git('ssh://nestle-eur@nestle-eur.visualstudio.com:22/NIHS_NEMO2/_git/NIHSnmr')

Quick start
-----------

    library("NIHSnmr")

### Play with a dataset

Load a dataset and interpolate it:

    nmrdata <- nmr_read_samples_dir("/dir/DUND-100713-Constitutive Thinness-plasma/",
                                    pulse_sequence = "NOESY") %>% nmr_interpolate()

View its metadata

    # Think of nmrdata as a list if you don't care about nmr_dataset objects
    View(nmrdata$metadata)

Get the matrix with the spectra:

    samples_matrix <- nmrdata$data_1r

Plot the spectra

    plot(nmrdata)

### Submit data to irods (at NIHS)

    sample_names <- c("dataset/10/", "dataset/20/", "dataset/30/") # samples to upload, you can as well use nm_read_samples_dir
    nmrdata <- nmr_read_samples(sample_names, metadata_only = TRUE)
    irods_metadata <- nmr_get_irods_meta(nmrdata)
    # Check irods_metadata and fill all the missing values
    View(irods_metadata)

    # Prepare zip files:
    irods_metadata <- nmr_prepare_zip_files(meta_irods = irods_metadata,
                                            workdir = 'my_dataset_zip_files')
    # Finally push the data to the irods directory
    nmr_push_to_irods(irods_metadata, "/NIHSData/DUND-XXXXXX/study/.../Metabolomics")

### Read samples from irods (at NIHS)

    # Search samples in irods
    irods_results <- nmr_irods_search(project_NPDI_ID = "DUND-100713",
                       study_nickname = "st_etienne",
                       assay_ID = "Plasma NMR metabolites",
                       assay_technique = "Nuclear Magnetic Resonance-1D-1H-NOESY")
    # Check that the results are the samples of interest
    View(irods_results)

    # Load the resulting irods samples:
    plasmadata <- nmr_read_samples_irods(irods_results)

See `?NIHSnmr` for further examples and links to documentation.

Graphical User Interface
------------------------

To open a graphical user interface you just need to run:

    library("NIHSnmr")
    NIHSnmr_interactive()

Contributing
------------

### At IBEC

To start contributing to this project you can use RStudio:

    File / New Project / Version Control / Git
    Repository URL: `git@gitnose.ibec.local:compmetabolomics/NIHSnmr.git`

### At NIHS

To start contributing to this project you can use RStudio:

    File / New Project / Version Control / Git
    Repository URL (choose):
     - ssh://nestle-eur@nestle-eur.visualstudio.com:22/NIHS_NEMO2/_git/NIHSnmr
     - https://nestle-eur.visualstudio.com/NIHS_NEMO2/_git/NIHSnmr
