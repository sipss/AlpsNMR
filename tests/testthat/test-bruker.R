## read_bruker_param / convert_field -------------------------------------------

test_that("read_bruker_param parses ParVal, Val continuations, Stamp, ParEmpty and Empty lines", {
    tmp <- withr::local_tempfile()
    writeLines(
        c(
            "##$TD= 65536",
            "##$SW_h= (0..1)",
            "1234.5",
            "##$BYTORDA= 0",
            "##TITLE= some title text",
            "$$ some stamp info",
            "##$EMPTYPAR=",
            "",
            "##$DTYPA= 2"
        ),
        tmp
    )

    result <- read_bruker_param(tmp)

    expect_equal(result$TD, 65536)
    # ParVec ("(0..1)") declares the field, and the following "Val" line's
    # value is attached to it and converted to numeric:
    expect_equal(result$SW_h, 1234.5)
    expect_equal(result$BYTORDA, 0)
    # A value that can't be parsed as numeric is kept as a string:
    expect_equal(result$TITLE, "some title text")
    expect_equal(result$Stamp, " some stamp info")
    # An empty parameter ("##$EMPTYPAR=" with nothing after it) is dropped:
    expect_false("EMPTYPAR" %in% names(result))
    expect_equal(result$DTYPA, 2)
})

test_that("read_bruker_param parses a vector value declared and given inline on the same line (ParVecVal)", {
    # e.g. "##$MYVEC= (0..2) 1 2 3": range declaration and values on one line,
    # as opposed to the more common case of a bare "(0..2)" declaration line
    # followed by the values on a separate line.
    tmp <- withr::local_tempfile()
    writeLines("##$MYVEC= (0..2) 1 2 3", tmp)

    result <- read_bruker_param(tmp)

    expect_equal(result$MYVEC, c(1, 2, 3))
})

test_that("read_bruker_param strips angle brackets from string values", {
    tmp <- withr::local_tempfile()
    writeLines("##$STR1= <hello world>", tmp)

    result <- read_bruker_param(tmp)

    expect_equal(result$STR1, "hello world")
})

test_that("convert_field converts numeric-looking strings to numeric vectors", {
    expect_equal(convert_field("123.5"), 123.5)
    expect_equal(convert_field("1 2 3"), c(1, 2, 3))
    expect_equal(convert_field("1e-3"), 0.001)
})

test_that("convert_field leaves non-numeric strings as strings", {
    expect_equal(convert_field("hello world"), "hello world")
    expect_equal(convert_field("1+e3"), "1+e3")
})

test_that("convert_field strips a single pair of surrounding angle brackets", {
    expect_equal(convert_field("<hello>"), "hello")
})

test_that("convert_field trims a whitespace-only string to an empty string", {
    expect_equal(convert_field(""), "")
    expect_equal(convert_field("  "), "")
})

## read_procs_file / read_acqus_file --------------------------------------------

test_that("read_procs_file reads only the procs files that exist by default", {
    tmp_dir <- withr::local_tempdir()
    writeLines("##$SI= 65536", file.path(tmp_dir, "procs"))
    writeLines("##$SI= 100", file.path(tmp_dir, "proc2s"))

    result <- read_procs_file(tmp_dir)

    expect_equal(names(result), c("procs", "proc2s"))
    expect_equal(result$procs$SI, 65536)
    expect_equal(result$proc2s$SI, 100)
})

test_that("read_procs_file honours an explicit procs_files argument", {
    tmp_dir <- withr::local_tempdir()
    writeLines("##$SI= 65536", file.path(tmp_dir, "procs"))
    writeLines("##$SI= 100", file.path(tmp_dir, "proc2s"))

    result <- read_procs_file(tmp_dir, procs_files = "procs")

    expect_equal(names(result), "procs")
})

test_that("read_acqus_file reads only the acqus files that exist by default", {
    tmp_dir <- withr::local_tempdir()
    writeLines("##$TD= 65536", file.path(tmp_dir, "acqus"))

    result <- read_acqus_file(tmp_dir)

    expect_equal(names(result), "acqus")
    expect_equal(result$acqus$TD, 65536)
})

## read_levels -------------------------------------------------------------------

test_that("read_levels applies each LEVSIGN case from a clevels text file", {
    make_clevels_dir <- function(levsign) {
        tmp_dir <- tempfile()
        dir.create(tmp_dir)
        writeLines(
            c(paste("##$LEVSIGN=", levsign), "##$MAXLEV= 3", "##$LEVELS= (0..5)", "1 -2 3 -4 5 -6"),
            file.path(tmp_dir, "clevels")
        )
        tmp_dir
    }

    expect_equal(read_levels(make_clevels_dir(1))$levels, c(1, 3, 5))
    expect_equal(read_levels(make_clevels_dir(2))$levels, c(-2, -4, -6))
    expect_equal(read_levels(make_clevels_dir(3))$levels, c(-2, -4, -6))
    expect_error(read_levels(make_clevels_dir(9)), "Unexpected clevels case")
})

test_that("read_levels duplicates a single remaining level into a length-2 vector", {
    tmp_dir <- withr::local_tempdir()
    writeLines(
        c("##$LEVSIGN= 1", "##$MAXLEV= 3", "##$LEVELS= (0..1)", "5 -6"),
        file.path(tmp_dir, "clevels")
    )

    result <- read_levels(tmp_dir)

    expect_equal(result$levels, c(5, 5))
})

test_that("read_levels returns NULL levels when neither level nor clevels exists", {
    tmp_dir <- withr::local_tempdir()
    expect_equal(read_levels(tmp_dir), list(levels = NULL))
})

test_that("read_levels reads the old binary level file given endian and NC_proc", {
    tmp_dir <- withr::local_tempdir()
    con <- file(file.path(tmp_dir, "level"), "wb")
    writeBin(as.integer(c(1, 3, 100, 200, 300)), con, size = 4, endian = "little")
    close(con)

    result <- read_levels(tmp_dir, endian = "little", NC_proc = 0)

    expect_equal(result$levels, c(100, 200, 300))
})

test_that("read_levels requires endian and NC_proc for the old binary level file", {
    tmp_dir <- withr::local_tempdir()
    con <- file(file.path(tmp_dir, "level"), "wb")
    writeBin(as.integer(c(1, 3, 100, 200, 300)), con, size = 4, endian = "little")
    close(con)

    expect_error(read_levels(tmp_dir), "endian information and NC_proc")
})

test_that("read_levels rejects a truncated binary level file", {
    tmp_dir <- withr::local_tempdir()
    con <- file(file.path(tmp_dir, "level"), "wb")
    writeBin(as.integer(c(1, 2)), con, size = 4, endian = "little")
    close(con)

    expect_error(
        read_levels(tmp_dir, endian = "little", NC_proc = 0),
        "not enough data"
    )
})

## guess_shape_and_submatrix_shape -----------------------------------------------

test_that("guess_shape_and_submatrix_shape returns all-NULL when procs info is missing or incomplete", {
    expect_equal(
        guess_shape_and_submatrix_shape(list()),
        list(dimension = NULL, shape = NULL, submatrix_shape = NULL)
    )
    expect_equal(
        guess_shape_and_submatrix_shape(list(procs = list(SI = 100))), # no XDIM
        list(dimension = NULL, shape = NULL, submatrix_shape = NULL)
    )
})

test_that("guess_shape_and_submatrix_shape infers 1D/2D/3D/4D shapes from the available proc*s entries", {
    base <- list(procs = list(SI = 100, XDIM = 32))
    expect_equal(guess_shape_and_submatrix_shape(base)[c("dimension", "shape", "submatrix_shape")],
        list(dimension = 1, shape = 100, submatrix_shape = 32))

    two_d <- c(base, list(proc2s = list(SI = 50, XDIM = 16)))
    expect_equal(guess_shape_and_submatrix_shape(two_d)[c("dimension", "shape", "submatrix_shape")],
        list(dimension = 2, shape = c(100, 50), submatrix_shape = c(32, 16)))

    three_d <- c(two_d, list(proc3s = list(SI = 20, XDIM = 8)))
    expect_equal(guess_shape_and_submatrix_shape(three_d)[c("dimension", "shape", "submatrix_shape")],
        list(dimension = 3, shape = c(100, 50, 20), submatrix_shape = c(32, 16, 8)))

    four_d <- c(three_d, list(proc4s = list(SI = 10, XDIM = 4)))
    expect_equal(guess_shape_and_submatrix_shape(four_d)[c("dimension", "shape", "submatrix_shape")],
        list(dimension = 4, shape = c(100, 50, 20, 10), submatrix_shape = c(32, 16, 8, 4)))
})

## parse_title_file ---------------------------------------------------------------

test_that("parse_title_file parses 'Name Value' pairs and trims trailing spaces/semicolons", {
    result <- parse_title_file(c("Sample John Doe ;", "Date 2024-01-01"))

    expect_equal(result$Sample, "John Doe")
    expect_equal(result$Date, "2024-01-01")
})

test_that("parse_title_file falls back to V1, V2... names when any line has no separate name/value", {
    result <- parse_title_file(c("Alanine", "Glucose"))

    expect_equal(result, list(V1 = "Alanine", V2 = "Glucose"))
})

## infer_dim_pulse_nuclei -----------------------------------------------------------

make_acqus_list <- function(exp, nuc1 = "1H", acqu2s_nuc1 = NULL, extra_nuclei = list()) {
    acqus <- c(list(EXP = exp, NUC1 = nuc1), extra_nuclei)
    for (n in setdiff(paste0("NUC", 2:8), names(extra_nuclei))) {
        acqus[[n]] <- "off"
    }
    out <- list(acqus = acqus)
    if (!is.null(acqu2s_nuc1)) {
        out$acqu2s <- list(NUC1 = acqu2s_nuc1)
    }
    out
}

test_that("infer_dim_pulse_nuclei recognizes each known pulse sequence and reports its nuclei", {
    expect_equal(infer_dim_pulse_nuclei(make_acqus_list("noesygppr1d"))$pulse_sequence, "NOESY")
    expect_equal(infer_dim_pulse_nuclei(make_acqus_list("cpmgpr1d"))$pulse_sequence, "CPMG")
    expect_equal(infer_dim_pulse_nuclei(make_acqus_list("diffSTE"))$pulse_sequence, "DIFFUSION")
    expect_equal(infer_dim_pulse_nuclei(make_acqus_list("jresgpprqf"))$pulse_sequence, "JRES")
    expect_equal(infer_dim_pulse_nuclei(make_acqus_list("zg30_PROTON"))$pulse_sequence, "PROTON")

    cosy <- infer_dim_pulse_nuclei(make_acqus_list("cosygpqf", acqu2s_nuc1 = "13C"))
    expect_equal(cosy$pulse_sequence, "COSY")
    expect_equal(cosy$nuclei, "1H-13C")
    expect_equal(cosy$dimension, 2) # acqus + acqu2s

    tocsy <- infer_dim_pulse_nuclei(make_acqus_list("mlevphpr", acqu2s_nuc1 = "1H"))
    expect_equal(tocsy$pulse_sequence, "TOCSY")
})

test_that("infer_dim_pulse_nuclei joins the non-'off' nuclei for HSQC/HMBC", {
    hsqc_acqus <- make_acqus_list("hsqcetgpsi", extra_nuclei = list(NUC2 = "13C"))
    hsqc <- infer_dim_pulse_nuclei(hsqc_acqus)
    expect_equal(hsqc$pulse_sequence, "HSQC")
    expect_equal(hsqc$nuclei, "1H-13C")

    hmbc_acqus <- make_acqus_list("hmbcgplpndqf", extra_nuclei = list(NUC2 = "13C"))
    hmbc <- infer_dim_pulse_nuclei(hmbc_acqus)
    expect_equal(hmbc$pulse_sequence, "HMBC")
    expect_equal(hmbc$nuclei, "1H-13C")
})

test_that("infer_dim_pulse_nuclei leaves pulse_sequence/nuclei as NA for an unrecognized experiment", {
    result <- infer_dim_pulse_nuclei(make_acqus_list("some-future-experiment"))

    expect_true(is.na(result$pulse_sequence))
    expect_true(is.na(result$nuclei))
    expect_equal(result$dimension, 1)
})

test_that("infer_dim_pulse_nuclei errors when the EXP field is missing", {
    expect_error(
        infer_dim_pulse_nuclei(list(acqus = list())),
        "EXP field is missing"
    )
})

## read_orig_file -----------------------------------------------------------------

test_that("read_orig_file returns NULL when there is no orig file", {
    tmp_dir <- withr::local_tempdir()
    expect_null(read_orig_file(tmp_dir))
})

test_that("read_orig_file parses 'name value' lines, joining multi-word values", {
    tmp_dir <- withr::local_tempdir()
    writeLines(
        c("Owner John Q Public", "Site my-site"),
        file.path(tmp_dir, "orig")
    )

    result <- read_orig_file(tmp_dir)

    expect_equal(result$Owner, "John Q Public")
    expect_equal(result$Site, "my-site")
})

## bruker_merge_meta_pdata ---------------------------------------------------------

test_that("bruker_merge_meta_pdata concatenates metadata and processed-data lists", {
    result <- bruker_merge_meta_pdata(list(a = 1), list(b = 2))
    expect_equal(result, list(a = 1, b = 2))
})
