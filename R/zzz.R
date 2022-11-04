# Hervé Pagès wrote:
#
# This is a race condition in 'R CMD check' that affects some packages on
# Windows. What 'R CMD check' is trying to do here is open the
# <package>-Ex.Rout file but that file is still being written to by another
# process. Windows (unlike Linux or macOS) doesn't like that: it won't let
# a process open a file if the file is being hold by another process.
# 
# More precisely here is what is happening here:
# 
# - One of the checks that 'R CMD check' performs is to run all the
# examples in a package, and to capture the output produced by the
# examples in a file (the <package>-Ex.Rout in your case).
# 
# - Once it's done running the examples, 'R CMD check' wants to know if
# the examples produced errors. It does this by parsing the
# <package>-Ex.Rout file. That's when it tries to open the file in read-only
# mode.
# 
# - But sometimes, when 'R CMD check' tries to open <package>-Ex.Rout in
# read-only mode, other processes are still holding on the <package>-Ex.Rout
# file. This typically happens when the examples in the package spawn
# subprocesses. For example these spawned processes can be the workers
# spawned by BiocParallel or by other parallel evaluation mechanism (e.g.
# doParallel/foreach). Note that there can be a very small lag between the
# moment the main worker and the spawned workers are done. If one of the
# spawned workers takes a fraction of an extra second to finish writing
# some output to <package>-Ex.Rout, then 'R CMD check' won't be able to open
# the file.
# 
# One could argue that this is an issue with 'R CMD check'. Maybe on
# Windows it should wait a couple of seconds before trying to open
# <package>-Ex.Rout, I don't know.
# 
# Maybe there's something that could be done in BiocParallel too. I don't
# know how hard that would be, or if that's even feasible, but maybe
# there's a way that the main worker could wait and make sure that all the
# spawn processes are dead before returning. I guess that kind of feature
# would only make sense for some parallel back-ends.
# 
# Does some of the examples in <package> use parallelization? If so maybe
# you can try to add a call to Sys.sleep(2) at the end of those examples
# and see if that helps.
# 
# In any case there's nothing that can really be done on our Windows
# builders to prevent this, sorry.

#' @rdname zzz
#' @name zzz
#' @title zzz
#' @keywords internal
#' @examples
#' # Workaround a bug in R CMD check
#' # Sys.sleep(2) # This worked, but this should be even more reliable:
#' BiocParallel::bpstop()
NULL
