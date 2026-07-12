## batman_get_full_filename ---------------------------------------------------

test_that("batman_get_full_filename creates the directory if missing and returns the joined path", {
    tmp <- file.path(withr::local_tempdir(), "does-not-exist-yet")
    expect_false(dir.exists(tmp))

    full_filename <- batman_get_full_filename(tmp, "foo.txt")

    expect_true(dir.exists(tmp))
    expect_equal(full_filename, file.path(tmp, "foo.txt"))
})

test_that("batman_get_full_filename errors if the target file already exists", {
    tmp <- withr::local_tempdir()
    file.create(file.path(tmp, "foo.txt"))

    expect_error(
        batman_get_full_filename(tmp, "foo.txt"),
        "already exists"
    )
})

## nmr_batman_options ---------------------------------------------------------

test_that("nmr_batman_options builds a batman_options object with the given ppmRange", {
    ppm_range <- matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE)
    bopts <- nmr_batman_options(ppmRange = ppm_range, specNo = "2")

    expect_s3_class(bopts, "batman_options")
    expect_equal(bopts$ppmRange, ppm_range)
    expect_equal(bopts$specNo, "2")
    # A couple of the numeric defaults, to catch accidental reordering of args:
    expect_equal(bopts$paraProc, 4L)
    expect_equal(bopts$nItBurnin, 200L)
})

test_that("nmr_batman_options rejects a ppmRange without exactly two columns", {
    expect_error(
        nmr_batman_options(ppmRange = matrix(1:3, ncol = 1)),
        "ncol"
    )
})

test_that("nmr_batman_options rejects a ppmRange containing NA", {
    expect_error(
        nmr_batman_options(ppmRange = matrix(c(1, NA, 2, 3), ncol = 2)),
        "anyNA"
    )
})

## nmr_batman_write_options ----------------------------------------------------

test_that("nmr_batman_write_options writes a options file and returns bopts unchanged", {
    tmp <- withr::local_tempdir()
    bopts <- nmr_batman_options(ppmRange = matrix(c(1, 2), ncol = 2), specNo = "1")

    result <- nmr_batman_write_options(bopts, batman_dir = tmp, filename = "opts.txt")

    expect_identical(result, bopts)
    written <- file.path(tmp, "opts.txt")
    expect_true(file.exists(written))
    lines <- readLines(written)
    expect_true(any(grepl("ppmRange - ppm ranges for analysis: \\(1,2\\)", lines)))
    expect_true(any(grepl("specNo - .*: 1$", lines)))
})

test_that("nmr_batman_write_options refuses to overwrite an existing options file", {
    tmp <- withr::local_tempdir()
    bopts <- nmr_batman_options(ppmRange = matrix(c(1, 2), ncol = 2))
    nmr_batman_write_options(bopts, batman_dir = tmp, filename = "opts.txt")

    expect_error(
        nmr_batman_write_options(bopts, batman_dir = tmp, filename = "opts.txt"),
        "already exists"
    )
})

## nmr_batman_export_dataset --------------------------------------------------

test_that("nmr_batman_export_dataset writes a ppm x NMRExperiment tab-separated table", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = c(1, 1.5, 2),
        data_1r = matrix(c(10, 20, 30, 40, 50, 60), nrow = 2, byrow = TRUE),
        metadata = list(external = data.frame(NMRExperiment = c("10", "20")))
    )
    tmp <- withr::local_tempdir()

    nmr_batman_export_dataset(dataset, batman_dir = tmp, filename = "NMRdata.txt")

    written <- file.path(tmp, "NMRdata.txt")
    expect_true(file.exists(written))
    table <- read.delim(written, check.names = FALSE)
    expect_equal(colnames(table), c("ppm", "NMRExperiment_10", "NMRExperiment_20"))
    expect_equal(table$ppm, c(1, 1.5, 2))
    expect_equal(table$NMRExperiment_10, c(10, 20, 30))
    expect_equal(table$NMRExperiment_20, c(40, 50, 60))
})

## nmr_batman_multi_data_user --------------------------------------------------

test_that("nmr_batman_multi_data_user adds default columns, sanitizes Metabolite and J_constant, and writes a csv", {
    multiplet_table <- data.frame(
        Metabolite = c("Ala,nine", "Glucose"),
        pos_in_ppm = c(1.5, 3.2),
        couple_code = c("0", "1"),
        J_constant = c(NA, 7.5),
        relative_intensity = c(1, 2)
    )
    tmp <- withr::local_tempdir()

    result <- nmr_batman_multi_data_user(multiplet_table, batman_dir = tmp, filename = "multi.csv")

    # Commas in metabolite names would break Batman's multiplet file format:
    expect_equal(result$Metabolite, c("Ala_nine", "Glucose"))
    # Batman does not accept NA for J_constant:
    expect_equal(result$J_constant, c(0, 7.5))
    # Default columns are added when missing:
    expect_equal(result$overwrite_pos, c(-50, -50))
    expect_equal(result$overwrite_truncation, c(-50, -50))
    expect_equal(result$Include_multiplet, c(1, 1))

    written <- read.csv(file.path(tmp, "multi.csv"))
    expect_equal(
        colnames(written),
        c(
            "Metabolite", "pos_in_ppm", "couple_code", "J_constant",
            "relative_intensity", "overwrite_pos", "overwrite_truncation",
            "Include_multiplet"
        )
    )
    expect_equal(written$Metabolite, c("Ala_nine", "Glucose"))
})

test_that("nmr_batman_multi_data_user preserves existing default-named columns instead of overwriting them", {
    multiplet_table <- data.frame(
        Metabolite = "Alanine",
        pos_in_ppm = 1.5,
        couple_code = "0",
        J_constant = 7.5,
        relative_intensity = 1,
        overwrite_pos = 99,
        overwrite_truncation = 99,
        Include_multiplet = 0
    )
    tmp <- withr::local_tempdir()

    result <- nmr_batman_multi_data_user(multiplet_table, batman_dir = tmp, filename = "multi.csv")

    expect_equal(result$overwrite_pos, 99)
    expect_equal(result$overwrite_truncation, 99)
    expect_equal(result$Include_multiplet, 0)
})

## nmr_batman_multi_data_user_hmdb --------------------------------------------

test_that("nmr_batman_multi_data_user_hmdb writes the bundled hmdb multiplet table", {
    tmp <- withr::local_tempdir()

    result <- nmr_batman_multi_data_user_hmdb(batman_dir = tmp, filename = "multi_hmdb.csv")

    expect_true(file.exists(file.path(tmp, "multi_hmdb.csv")))
    expect_true("Metabolite" %in% colnames(result))
    expect_gt(nrow(result), 0)
})

## nmr_batman_metabolites_list -------------------------------------------------

test_that("nmr_batman_metabolites_list writes unique metabolite names and returns the input unchanged", {
    tmp <- withr::local_tempdir()
    metabolite_names <- c("alanine", "glucose", "alanine")

    result <- nmr_batman_metabolites_list(metabolite_names, batman_dir = tmp, filename = "metsList.csv")

    # The return value is the original (non-deduplicated) input:
    expect_equal(result, metabolite_names)
    # But the written file is deduplicated:
    written <- readLines(file.path(tmp, "metsList.csv"))
    expect_equal(written, c("alanine", "glucose"))
})
