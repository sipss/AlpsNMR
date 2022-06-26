#' Download MTBLS242
#' 
#' Downloads part of the MTBLS242 dataset (few hundred megabytes)
#' 
#' NMR dataset from Gralka et al., 2015. DOI: 10.3945/ajcn.115.110536.
#'
#' @param dest_dir Directory where the dataset should be saved
#' @param force Logical. If `TRUE` we do not re-download files if they exist
#'
#' @return A data frame with the downloaded samples
#' @export
#'
#' @examples
#' \dontrun{
#' download_MTBLS242("./MTBLS242")
#' }
download_MTBLS242 <- function(dest_dir = "MTBLS242", force = FALSE) {
    require_pkgs(pkg = c("curl", "zip", "archive"))
    url <- "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS242/"
    
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Download metadata file s_mtbls242.txt
    meta_file <- "s_mtbls242.txt"
    meta_url <- file.path(url, meta_file)
    
    # meta_dst <- file.path(dest_dir, meta_file)
    # utils::download.file(meta_url, method = "auto", destfile = meta_dst, mode = "wb")
    annotations_destfile <- file.path(dest_dir, "sample_annotations.tsv")
    if (!file.exists(annotations_destfile) || force) {
        rlang::inform("Downloading sample annotations, keeping only preop and 3 months after surgery timepoints...")
        sample_annot <- tibble::as_tibble(
            utils::read.table(
                meta_url,
                sep = "\t",
                header = TRUE,
                check.names = FALSE
            ),
            .name_repair = "minimal"
        )
        
        # Keep sample name and time point:
        sample_annot <- dplyr::select(
            sample_annot,
            c("NMRExperiment" = "Sample Name", "TimePoint" = "Factor Value[time point]")
        )
        sample_annot <- dplyr::filter(sample_annot, .data$TimePoint %in% c("preop", "3 months after surgery"))
        sample_annot$NMRExperiment <- gsub(pattern = "-", replacement = "_", sample_annot$NMRExperiment, fixed = TRUE)
        
        sample_annot <- tidyr::separate(
            sample_annot,
            col = "NMRExperiment",
            into = c("timepoint", "SampleID", "S"),
            sep = "_",
            remove = FALSE
        )
        sample_annot <- dplyr::select(sample_annot, -"S")
        sample_annot$filename <- paste0("Obs", sample_annot$timepoint, "_", sample_annot$SampleID, "s")
        
        # File Obs1_0256s.zip incorrectly contains Obs1_0010s. Remove that ID (in all timepoints)
        sample_annot <- dplyr::filter(sample_annot, .data$SampleID != "0256")
        
        # Keep samples matched in the two timepoints under study:
        sample_annot <- dplyr::group_by(sample_annot, .data$SampleID)
        sample_annot <- dplyr::filter(sample_annot, dplyr::n() == 2)
        sample_annot <- dplyr::ungroup(sample_annot)
        sample_annot$NMRExperiment <- sample_annot$filename
        sample_annot <- dplyr::select(sample_annot, -"timepoint", -"filename")
        
        utils::write.table(sample_annot, file = annotations_destfile, sep = "\t", row.names = FALSE)
    } else {
        rlang::inform(glue("Annotations were previously downloaded. Loading {annotations_destfile}"))
        sample_annot <- readr::read_tsv(annotations_destfile)
    }
    rlang::inform("Downloading samples and keeping CPMG spectra, please wait...")
    pb <- progress_bar_new(
        name = "Preparing samples",
        total = length(sample_annot$NMRExperiment)
    )
    dst_rootdir <- file.path(dest_dir, "samples")
    dir.create(dst_rootdir, recursive = TRUE, showWarnings = FALSE)
    purrr::walk(
        sample_annot$NMRExperiment,
        function(filename_base, url, dst_rootdir, keep_only_CPMG_1r) {
            progress_bar_update(pb)
            filename <- paste0(filename_base, ".zip")
            src_url <- file.path(url, filename)
            final_dst_file <- file.path(dst_rootdir, filename)
            intermediate_dst_file <- file.path(dst_rootdir, paste0(filename, "intermediate.zip"))
            if (file.exists(final_dst_file) && !force) {
                rlang::inform(glue("Skipping download of {filename_base} since it was previously downloaded."))
                return()
            }
            curl::curl_download(url = src_url, destfile = intermediate_dst_file)
            if (keep_only_CPMG_1r) {
                # zip file management in R is as of mid 2022 a pain:
                # Base R provides tools that are not portable (system dependent)
                # The `zip` package segfaults (to me and to many github issues)
                # The `archive` package does not provide an API to extract selectively folders from a zip file.
                # We'll see if the situation improves in some years, so we can clean this mess:
                # For now, we unzip with zip::unzip() because it is portable, does not segfault and allows
                # us to extract a subdir.
                # We will zip with archive because it is portable and does not segfault
                filenames_in_zip <- zip::zip_list(intermediate_dst_file)[["filename"]]
                prefix_to_keep <- file.path(filename_base, "3", "") # subdirectory 3/ contains the CPMG measurement
                filenames_in_zip <- filenames_in_zip[startsWith(filenames_in_zip, prefix_to_keep)]
                zip::unzip(zipfile = intermediate_dst_file, exdir = dst_rootdir, files = filenames_in_zip)
                unlink(intermediate_dst_file)
                cwd <- getwd()
                setwd(dst_rootdir)
                tempname <- paste0(filename_base, "-rm")
                file.rename(filename_base, tempname)
                file.rename(file.path(tempname, "3"), filename_base)
                unlink(tempname, recursive = TRUE)
                unlink(file.path(filename_base, "fid"))
                unlink(file.path(filename_base, "pdata", "1", "1i"))
                # pick your poison:
                # a) utils::zip() is annoyingly verbose and possibly not portable
                #utils::zip(zipfile = filename, files = filename_base)
                #setwd(cwd)
                # b) zip::zip() segfaults
                #zip::zip(zipfile = final_dst_file, files = filename_base, root = dst_rootdir)
                # c) archive (yet another package)
                archive::archive_write_dir(archive = filename, dir = filename_base)
                setwd(cwd)
                # And once you have the zip file, remove the directory:
                unlink(file.path(dst_rootdir, filename_base), recursive = TRUE)
            }
        },
        url = url,
        dst_rootdir = dst_rootdir,
        keep_only_CPMG_1r = TRUE
    )
    progress_bar_end(pb)
    sample_annot
}
