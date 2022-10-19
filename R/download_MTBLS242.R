#' Download MTBLS242
#'
#' Downloads the [MTBLS242](https://www.ebi.ac.uk/metabolights/MTBLS242/protocols)
#' dataset from Gralka et al., 2015. DOI: \doi{10.3945/ajcn.115.110536}.
#' 
#' Besides the destination directory, this function
#' includes three logical parameters to limit the amount of downloaded/saved data.
#' To run the tutorial workflow:
#' - only the "preop" and "three months" timepoints are used,
#' - only subjects measured in *both* preop and three months time points are used
#' - only the CPMG samples are used.
#' 
#' If you want to run the tutorial, you can set those filters to `TRUE`. Then, roughly
#' 800MB will be downloaded, and 77MB of disk space will be used, since for each
#' downloaded sample we remove all the data but the CPMG.
#' 
#' If you set those filters to `FALSE`, roughly 1.8GB of data will be
#' downloaded (since we have more timepoints to download) and 1.8GB
#' of disk space will be used.
#' 
#' Note that we have experienced some sporadic timeouts from Metabolights, 
#' when downloading the dataset. If you get those timeouts simply re-run the
#' download function and it will restart from where it stopped.
#' 
#' Note as well, that we observed several files to have incorrect data:
#' - Obs4_0346s.zip is not present in the FTP server
#' - Obs0_0110s.zip and Obs1_0256s.zip incorrectly contain sample Obs1_0010s
#' 
#' This function removes all three samples from the samples annotations and
#' doesn't download their data.
#' 
#' 
#' @param dest_dir Directory where the dataset should be saved
#' @param force Logical. If `TRUE` we do not re-download files if they exist. The function does not check whether cached versions were
#' downloaded with different `keep_only_*` arguments, so please use `force = TRUE` if you change the `keep_only_*` settings.
#' @param keep_only_CPMG_1r If `TRUE`, remove all other data beyond the CPMG real spectrum, which is enough for the tutorial
#' @param keep_only_preop_and_3months If `TRUE`, keep only the preoperatory and the "three months after surgery" time points, enough for the tutorial
#' @param keep_only_complete_time_points If `TRUE`, remove samples that do not appear on all timepoints. Useful for the tutorial.
#'
#' @return Invisibly, the annotations. See the example for how to download the
#'  annotations and create a dataset from the downloaded files.
#' @export
#'
#' @examples
#' \dontrun{
#'   download_MTBLS242("./MTBLS242")
#'   annot <- readr::read_tsv(annotations_destfile)
#'   
#'   dataset <- nmr_read_samples(annot$filename)
#'   dataset <- nmr_meta_add(dataset, annot)
#'   dataset
#' }
download_MTBLS242 <- function(
        dest_dir = "MTBLS242", force = FALSE,
        keep_only_CPMG_1r = TRUE,
        keep_only_preop_and_3months = TRUE,
        keep_only_complete_time_points = TRUE
    ) {
    require_pkgs(pkg = c("curl", "zip", "archive"))
    url <- "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS242/"

    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

    # Download metadata file s_mtbls242.txt
    meta_file <- "s_mtbls242.txt"
    meta_url <- file.path(url, meta_file)

    # meta_dst <- file.path(dest_dir, meta_file)
    # utils::download.file(meta_url, method = "auto", destfile = meta_dst, mode = "wb")
    annotations_destfile <- file.path(dest_dir, "sample_annotations.tsv")
    annotations_orig_destfile <- file.path(dest_dir, "s_mtbls242.txt")
    dst_rootdir <- file.path(dest_dir, "samples")
    if (!file.exists(annotations_orig_destfile) || force) {
        rlang::inform(c("i" = "Downloading sample annotations..."))
        curl::curl_download(url = meta_url, destfile = annotations_orig_destfile)
    }
    if (!file.exists(annotations_destfile) || force) {
        sample_annot <- tibble::as_tibble(
            utils::read.table(
                annotations_orig_destfile,
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
        if (keep_only_preop_and_3months) {
            sample_annot <- dplyr::filter(sample_annot, .data$TimePoint %in% c("preop", "3 months after surgery"))
        }
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

        sample_annot$NMRExperiment <- sample_annot$filename
        sample_annot <- dplyr::select(sample_annot, -"timepoint", -"filename")
        
        # File Obs0_0110s.zip incorrectly contains Obs0_0010s. Remove that ID
        sample_annot <- dplyr::filter(sample_annot, .data$NMRExperiment != "Obs0_0110s")
        # File Obs1_0256s.zip incorrectly contains Obs1_0010s. Remove that ID
        sample_annot <- dplyr::filter(sample_annot, .data$NMRExperiment != "Obs1_0256s")
        # File Obs4_0346s.zip does not exist in the FTP server, remove that entry:
        sample_annot <- dplyr::filter(sample_annot, .data$NMRExperiment != "Obs4_0346s")

        
        # Keep samples matched in the two timepoints under study:
        num_timepoints <- length(unique(sample_annot$TimePoint))
        if (keep_only_complete_time_points) {
            sample_annot <- dplyr::group_by(sample_annot, .data$SampleID)
            sample_annot <- dplyr::filter(sample_annot, dplyr::n() == !!num_timepoints)
            sample_annot <- dplyr::ungroup(sample_annot)
        }
        
        # filename.zip!/path/to/sample/in/zip
        if (keep_only_CPMG_1r) {
            sample_annot$filename <- paste0(
                dst_rootdir, "/", sample_annot$NMRExperiment, ".zip",
                "!",
                "/", sample_annot$NMRExperiment
            )
        } else {
            sample_annot$filename <- paste0(
                dst_rootdir, "/", sample_annot$NMRExperiment, ".zip",
                "!",
                "/", sample_annot$NMRExperiment, "/3"
            )
            
        }
        utils::write.table(sample_annot, file = annotations_destfile, sep = "\t", row.names = FALSE)
    } else {
        rlang::inform(c("i" = glue("Annotations were previously saved. Loading {annotations_destfile}")))
        sample_annot <- utils::read.csv(annotations_destfile, header = TRUE, sep = "\t")
    }
    rlang::inform(c("i" = "Downloading samples, please wait..."))
    pb <- progress_bar_new(
        name = "Preparing samples",
        total = length(sample_annot$NMRExperiment)
    )
    dir.create(dst_rootdir, recursive = TRUE, showWarnings = FALSE)
    report_skipped_downloads <- FALSE
    purrr::walk(
        sample_annot$NMRExperiment,
        function(filename_base, url, dst_rootdir, keep_only_CPMG_1r) {
            progress_bar_update(pb)
            filename <- paste0(filename_base, ".zip")
            src_url <- file.path(url, filename)
            final_dst_file <- file.path(dst_rootdir, filename)
            intermediate_dst_file <- file.path(dst_rootdir, paste0(filename, "intermediate.zip"))
            if (file.exists(final_dst_file) && !force) {
                if (!report_skipped_downloads) {
                    rlang::inform(c("i" = "Skipping re-download of previously downloaded samples."))
                    report_skipped_downloads <<- TRUE
                }
                return()
            }
            tryCatch(
                {
                    curl::curl_download(url = src_url, destfile = intermediate_dst_file)
                },
                error = function(e) {
                    msg <- conditionMessage(e)
                    rlang::abort(
                        message = c(
                            "Sample failed to download",
                            "i" = glue("Sample: {filename_base}"),
                            "i" = glue("URL: {src_url}"),
                            "i" = glue("Original error message: {msg}")
                        ),
                        parent = e
                    )
                }
            )
            if (!keep_only_CPMG_1r) {
                file.rename(intermediate_dst_file, final_dst_file)
            } else {
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
                # utils::zip(zipfile = filename, files = filename_base)
                # setwd(cwd)
                # b) zip::zip() segfaults
                # zip::zip(zipfile = final_dst_file, files = filename_base, root = dst_rootdir)
                # c) archive (yet another package)
                # I have to put everything in a folder:
                tmpname <- paste0(filename_base, "-base")
                dir.create(tmpname)
                file.rename(filename_base, file.path(tmpname, filename_base))
                file.rename(tmpname, filename_base)
                archive::archive_write_dir(archive = filename, dir = filename_base)
                setwd(cwd)
                # And once you have the zip file, remove the directory:
                unlink(file.path(dst_rootdir, filename_base), recursive = TRUE)
            }
        },
        url = url,
        dst_rootdir = dst_rootdir,
        keep_only_CPMG_1r = keep_only_CPMG_1r
    )
    progress_bar_end(pb)
    invisible(sample_annot)
}
