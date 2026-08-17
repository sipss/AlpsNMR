cat("== TBB/Rfast diagnostic ==\n")
cat("R.version.string:", R.version.string, "\n")
cat("HTTPUserAgent:", getOption("HTTPUserAgent"), "\n")
cat("options(repos):\n"); print(getOption("repos"))
cat("Sys.getenv('RSPM'):", Sys.getenv("RSPM"), "\n")

ip <- installed.packages()
for (pkg in c("Rfast", "RcppParallel")) {
  if (pkg %in% rownames(ip)) {
    cat(pkg, "version:", ip[pkg, "Version"], "Built:", ip[pkg, "Built"],
        "LibPath:", ip[pkg, "LibPath"], "\n")
  } else cat(pkg, "NOT installed\n")
}

# Every libtbb on the filesystem, and whether it still exports the classic
# tbb::task typeinfo Rfast.so was (presumably) linked against.
libs <- system("find / -xdev -name 'libtbb*.so*' 2>/dev/null", intern = TRUE)
cat("libtbb libraries found:\n"); print(libs)
for (lib in libs) {
  has_task <- system(sprintf("nm -D '%s' 2>/dev/null | grep -q _ZTIN3tbb4taskE && echo yes || echo no", lib), intern = TRUE)
  cat(" ", lib, "-> exports tbb::task:", has_task, "\n")
}

if ("Rfast" %in% rownames(ip)) {
  so <- file.path(find.package("Rfast"), "libs", "Rfast.so")
  cat("ldd", so, ":\n"); system(paste("ldd", shQuote(so)))

  # Does a bare `library(Rfast)` already fail, outside any parallel context?
  r1 <- system2("Rscript", c("-e", shQuote("library(Rfast)")), stdout = TRUE, stderr = TRUE)
  cat("plain library(Rfast) in fresh subprocess:\n", paste(r1, collapse = "\n"), "\n")

  # Does the actual bplapply pattern fail here, right now?
  ok <- tryCatch({
    library(BiocParallel)
    print(bplapply(1:2, function(i) { library(Rfast); Rfast::colmeans(matrix(rnorm(100), ncol = 10), parallel = TRUE) }))
    TRUE
  }, error = function(e) { cat("bplapply FAILED:", conditionMessage(e), "\n"); FALSE })
  cat("bplapply reproduction succeeded:", ok, "\n")

  # Is the CRAN mirror actually serving a binary, and for what platform?
  ver <- as.character(packageVersion("Rfast"))
  url <- paste0(Sys.getenv("RSPM"), "/src/contrib/Rfast_", ver, ".tar.gz")
  system(sprintf("curl -sI -A %s '%s' | grep -iE 'x-package-type|x-package-binary-tag|^HTTP'",
                  shQuote(getOption("HTTPUserAgent")), url))
}
