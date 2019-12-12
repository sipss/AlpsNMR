#' Downloadind the MTBLS242 dataset
#'
#' @param to A directory
#'
#' @return A folder with demo samples
#' @export
#' @example 
#' \dontrun{
#' download_demo(to = "C:/Users/")
#'}
download_demo <- function(to = ".") {
  if (!dir.exists(to)) {
    dir.create(to, recursive = TRUE)
  }
  zip_fn <- file.path(to, "MTBLS242.zip")
  dropbox_link <- "https://dl.dropboxusercontent.com/s/0snivrsd7m82yey/MTBLS242.zip?dl=0"
  utils::download.file(dropbox_link, method = "auto", destfile = zip_fn, mode = "wb")
  
  # output <- file.path(to, "MTBLS242")
  # fs::dir_create(output)
  
  zip::unzip(zipfile = zip_fn, exdir = to)
  invisible(NULL)
}