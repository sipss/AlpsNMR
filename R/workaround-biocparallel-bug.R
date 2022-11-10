# This whole file is a workaround for:
# - https://github.com/Bioconductor/BiocParallel/pull/227
#   (merged released on Bioconductor 3.16, BiocParallel 1.31.14)
# And a workaround for:
# - https://github.com/Bioconductor/BiocParallel/pull/228
#   (as of 2022-11-04 it is under review)
#
# Please once that's merged and released, ensure you have
# BiocParallel (>= 1.xx.xx?)
# in the DESCRIPTION file, replace all usages of mymapply() with
# BiocParallel::bpmapply() and delete this file.
#
# Most of these functions are copies from the BiocParallel package

.getDotsForMapply <- function (...)
{
  ddd <- list(...)
  if (!length(ddd))
    return(list(list()))
  len <- vapply(ddd, length, integer(1L))
  if (!all(len == len[1L])) {
    max.len <- max(len)
    if (max.len && any(len == 0L))
      stop("zero-length and non-zero length inputs cannot be mixed")
    if (any(max.len%%len))
      warning("longer argument not a multiple of length of vector")
    ddd <- lapply(ddd, rep_len, length.out = max.len)
  }
  ddd
}

.simplify <- function (results, SIMPLIFY = FALSE)
{
  if (SIMPLIFY && length(results))
    results <- simplify2array(results)
  results
}

# see test_utilities.R:test_transposeArgsWithIterations() for all
# USE.NAMES corner cases
.transposeArgsWithIterations <- function(nestedList, USE.NAMES) {
  num_arguments <- length(nestedList)
  if (num_arguments == 0L) {
    return(list())
  }

  # nestedList[[1L]] has the values for the first argument in all iterations
  num_iterations <- length(nestedList[[1L]])

  # count the iterations, and name them if needed
  iterations <- seq_len(num_iterations)
  if (USE.NAMES) {
    first_arg <- nestedList[[1L]]
    if (is.character(first_arg) && is.null(names(first_arg))) {
      names(iterations) <- first_arg
    } else {
      names(iterations) <- names(first_arg)
    }
  }

  # argnames:
  argnames <- names(nestedList)

  # on iteration `i` we get the i-th element from each list.
  # note that .getDotsForMapply() has taken care already of ensuring that
  # nestedList elements are recycled properly
  lapply(iterations, function(i) {
    x <- lapply(nestedList, function(argi) {
      unname(argi[i])
    })
    names(x) <- argnames
    x
  })
}

.wrapMapply <- function(dots, .FUN, .MoreArgs) {
  .mapply(.FUN, dots, .MoreArgs)[[1L]]
}

mymapply <- function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE,
                     BPREDO = list(), BPPARAM = BiocParallel::bpparam(), BPOPTIONS = BiocParallel::bpoptions()) {
  ## re-package for lapply
  ddd <- .getDotsForMapply(...)
  if (!length(ddd))
    return(list())

  ddd <- .transposeArgsWithIterations(ddd, USE.NAMES)
  if (!length(ddd))
    return(ddd)

  FUN <- match.fun(FUN)

  res <- BiocParallel::bplapply(X=ddd, .wrapMapply, .FUN=FUN,
                  .MoreArgs=MoreArgs, BPREDO=BPREDO,
                  BPPARAM=BPPARAM, BPOPTIONS = BPOPTIONS)
  .simplify(res, SIMPLIFY)
}


