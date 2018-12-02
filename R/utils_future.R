# Reexport some functions from the future package for our confort

#' @importFrom future plan sequential multicore multiprocess multisession cluster remote
#' @export
future::plan

#' @export
future::sequential

#' @export
future::multicore

#' @export
future::multiprocess

#' @export
future::multisession

#' @export
future::cluster

#' @export
future::remote
