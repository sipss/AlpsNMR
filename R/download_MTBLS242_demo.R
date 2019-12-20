#' Download the MTBLS242 dataset
#' 
#' The function downloads the NMR dataset from Gralka et al., 2015. DOI: 10.3945/ajcn.115.110536.
#'
#' @param to A directory
#' 
#' @return A folder with demo samples
#' @export
#' @references \url{https://doi.org/10.3945/ajcn.115.110536}
#' 
#' @examples 
#' \dontrun{
#' library(AlpsNMR)
#' download_demo(to = "C:/Users/")
#' }
#' 
download_demo <- function(to = ".") {
  if (!dir.exists(to)) {
    dir.create(to, recursive = TRUE)
  }
  zip_fn <- file.path(to, "MTBLS242.zip")
  dropbox_link <- "https://dl.dropboxusercontent.com/s/0snivrsd7m82yey/MTBLS242.zip?dl=0"
  utils::download.file(dropbox_link, method = "auto", destfile = zip_fn, mode = "wb")
  zip::unzip(zipfile = zip_fn, exdir = to)
  invisible(NULL)
}