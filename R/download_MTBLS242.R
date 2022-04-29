#' Download MTBLS242
#' 
#' Downloads part of the MTBLS242 dataset (few hundred megabytes)
#' 
#' NMR dataset from Gralka et al., 2015. DOI: 10.3945/ajcn.115.110536.
#'
#' @param dest_dir Directory where the dataset should be saved
#'
#' @return A data frame with the downloaded samples
#' @export
#'
#' @examples
#' \dontrun{
#' download_MTBLS242("./MTBLS242")
#' }
download_MTBLS242 <- function(dest_dir = "MTBLS242") {
  require_pkgs(pkg = c("curl", "zip"))
  url <- "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS242/"
  
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Downloading sample annotations, keeping only preop and 3 months after surgery timepoints...")
  # Download metadata file s_mtbls242.txt
  meta_file <- "s_mtbls242.txt"
  meta_url <- file.path(url, meta_file)
  
  # meta_dst <- file.path(dest_dir, meta_file)
  # utils::download.file(meta_url, method = "auto", destfile = meta_dst, mode = "wb")
  
  sample_annot <- tibble::as_tibble(utils::read.table(meta_url, sep = "\t"))
  
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
  sample_annot <- dplyr::select(sample_annot, -"timepoint")
  
  utils::write.table(sample_annot, file = file.path(dest_dir, "sample_annotations.tsv"), sep = "\t", row.names = FALSE)
  
  message("Downloading samples and keeping CPMG spectra, please wait...")
  pb <- progress_bar_new(
    name = "Preparing samples",
    total = length(sample_annot$filename)
  )
  dst_rootdir <- file.path(dest_dir, "samples")
  dir.create(dst_rootdir, recursive = TRUE, showWarnings = FALSE)
  purrr::walk(
    sample_annot$filename,
    function(filename_base, url, dst_rootdir, keep_only_CPMG = FALSE, progress_bar = NULL) {
      progress_bar_update(pb = pb)
      filename <- paste0(filename_base, ".zip")
      src_url <- file.path(url, filename)
      dst_file <- file.path(dst_rootdir, filename)
      if (file.exists(dst_file)) {
        return()
      }
      curl::curl_download(url = src_url, destfile = dst_file)
      if (keep_only_CPMG) {
        filenames_in_zip <- zip::zip_list(dst_file)[["filename"]]
        prefix_to_keep <- file.path(filename_base, "3", "") # subdirectory 3/ contains the CPMG measurement
        filenames_in_zip <- filenames_in_zip[startsWith(filenames_in_zip, prefix_to_keep)]
        zip::unzip(zipfile = dst_file, exdir = dst_rootdir, files = filenames_in_zip)
        unlink(dst_file)
        zip::zip(zipfile = dst_file, files = filename_base, root = dst_rootdir)
        unlink(file.path(dst_rootdir, filename_base), recursive = TRUE)
      }
    },
    url = url,
    dst_rootdir = dst_rootdir,
    keep_only_CPMG = TRUE,
    progress_bar = pb
  )
  progress_bar_end(pb)
  sample_annot
}
