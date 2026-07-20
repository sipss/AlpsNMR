# Download MTBLS242

Downloads the
[MTBLS242](https://www.ebi.ac.uk/metabolights/MTBLS242/protocols)
dataset from Gralka et al., 2015. DOI:
[doi:10.3945/ajcn.115.110536](https://doi.org/10.3945/ajcn.115.110536) .

## Usage

``` r
download_MTBLS242(
  dest_dir = "MTBLS242",
  force = FALSE,
  keep_only_CPMG_1r = TRUE,
  timepoints = c("preop", "12 months after surgery"),
  keep_only_complete_time_points = TRUE
)
```

## Arguments

- dest_dir:

  Directory where the dataset should be saved. Every freshly downloaded
  file is verified against the canonical SHA-256 checksums MetaboLights
  publishes for MTBLS242. As a fallback for when that manifest cannot be
  fetched, the SHA-256 of every downloaded file is also pinned to
  `<dest_dir>/SHA256SUMS` the first time it is saved, and re-verified on
  every later call that reuses a cached file, so local corruption or
  tampering between calls is detected either way.

- force:

  Logical. If `TRUE` we do not re-download files if they exist. The
  function does not check whether cached versions were downloaded with
  different `keep_only_*`/`timepoints` arguments, so please use
  `force = TRUE` if you change those settings. `force = TRUE` also
  re-downloads and re-pins the checksum of every file, rather than
  verifying it against a previously pinned value. If you only *added*
  timepoints (kept every previously requested one), deleting
  `<dest_dir>/sample_annotations.tsv` before calling again is enough to
  pick up the new timepoints without re-downloading previously cached
  samples.

- keep_only_CPMG_1r:

  If `TRUE`, remove all other data beyond the CPMG real spectrum, which
  is enough for the tutorial

- timepoints:

  Either `NULL` to keep every timepoint, or a character vector with the
  `TimePoint` values to keep (a subset of `"preop"`,
  `"3 months after surgery"`, `"6 months after surgery"`,
  `"9 months after surgery"`, `"12 months after surgery"`). Defaults to
  `c("preop", "12 months after surgery")`, enough for the tutorial.

- keep_only_complete_time_points:

  If `TRUE`, remove samples that do not appear on all of the kept
  timepoints. Useful for the tutorial.

## Value

Invisibly, the annotations. See the example for how to download the
annotations and create a dataset from the downloaded files.

## Details

Besides the destination directory, this function includes parameters to
limit the amount of downloaded/saved data. To run the tutorial workflow
with a two-timepoint comparison:

- only the requested `timepoints` (e.g. "preop" and "12 months after
  surgery") are used,

- only subjects measured in *all* of the requested `timepoints` are
  used,

- only the CPMG samples are used.

The full study has five timepoints: `"preop"`,
`"3 months after surgery"`, `"6 months after surgery"`,
`"9 months after surgery"`, and `"12 months after surgery"`. Restricting
`timepoints` to two of them keeps the download small (roughly 800MB
downloaded, 77MB on disk with `keep_only_CPMG_1r = TRUE`). Passing
`timepoints = NULL` downloads every timepoint (roughly 1.8GB downloaded
and on disk).

Note that we have experienced some sporadic timeouts from Metabolights,
when downloading the dataset. If you get those timeouts simply re-run
the download function and it will restart from where it stopped.

Note as well, that we observed several files to have incorrect data:

- Obs4_0346s.zip is not present on the server

- Obs0_0110s.zip and Obs1_0256s.zip incorrectly contain sample
  Obs1_0010s

This function removes all three samples from the samples annotations and
doesn't download their data.

## Examples

``` r
if (FALSE) { # \dontrun{
  download_MTBLS242("./MTBLS242")
  annot <- readr::read_tsv(annotations_destfile)

  dataset <- nmr_read_samples(annot$filename)
  dataset <- nmr_meta_add(dataset, annot)
  dataset
} # }
```
