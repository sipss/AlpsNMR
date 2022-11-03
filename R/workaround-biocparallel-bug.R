# This whole file is a workaround for:
# - https://github.com/Bioconductor/BiocParallel/pull/227
# And a workaround for:
# - https://github.com/Bioconductor/BiocParallel/pull/228
#
# Here, for #228, if bpparam() is multicore or serial, memory is shared and
# we use the old approach. If memory is not shared we use the proposed approach.
# So we get the best of both worlds. I wonder what will be decided on #228
# and whether there will be a trade-off or not
#
# Please once that's merged and released, ensure you have
# BiocParallel (>= 1.xx.xx?)
# in the DESCRIPTION file, replace all usages of mymapply() with
# BiocParallel::bpmapply() and delete this file.
#
# Some of these functions are copies from the BiocParallel package
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

mybpmapply_not_share <- function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE,
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

## re-apply the names on, e.g., mapply(FUN, X) to the return value;
## see test_mrename for many edge cases
.mrename <- function(results, dots, USE.NAMES=TRUE) {
    ## dots: a list() containing one element for each ... argument
    ## passed to mapply
    if (USE.NAMES) {
      ## extract the first argument; if there are no arguments, then
      ## dots is (unnamed) list(0)
      if (length(dots))
        dots <- dots[[1L]]
      if (is.character(dots) && is.null(names(dots))) {
        names(results) <- dots
      } else {
        names(results) <- names(dots)
      }
    } else {
      results <- unname(results)
    }
    results
}

.wrap <- function(.i, .FUN, .ddd, .MoreArgs) {
  dots <- lapply(.ddd, `[`, .i)
  .mapply(.FUN, dots, .MoreArgs)[[1L]]
}

mybpmapply_share <- function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE,
                             BPREDO = list(), BPPARAM = BiocParallel::bpparam(), BPOPTIONS = BiocParallel::bpoptions()) {
  ## re-package for lapply
  ddd <- .getDotsForMapply(...)
  if (!length(ddd) || !length(ddd[[1L]]))
    return(.mrename(list(), ddd, USE.NAMES))

  FUN <- match.fun(FUN)

  res <- BiocParallel::bplapply(X=seq_along(ddd[[1L]]), .wrap, .FUN=FUN, .ddd=ddd,
                  .MoreArgs=MoreArgs, BPREDO=BPREDO,
                  BPPARAM=BPPARAM, BPOPTIONS = BPOPTIONS)
  .simplify(.mrename(res, ddd, USE.NAMES), SIMPLIFY)
}

workers_share_memory <- function() {
  tryCatch({
    bp_param <- BiocParallel::bpparam()
    methods::is(bp_param, "SerialParam") || methods::is(bp_param, "MulticoreParam")
  },
  error = function(e) TRUE
  )
}

mymapply <- function(...) {
  # If workers can share memory, it should be better to pass a copy of all arguments
  # to all workers.
  # Otherwise, work on transposing the arguments so each worker gets its own
  # slice to work.
  if (workers_share_memory()) {
    mybpmapply_share(...)
  } else {
    mybpmapply_not_share(...)
  }
}

