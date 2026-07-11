test_that("nmr_dataset_autophase works", {
    skip_if_not_installed("NMRphasing")
    lorentzian <- function(x, x0, gamma, A) {
        A * (1 / (pi * gamma)) * ((gamma^2) / ((x - x0)^2 + gamma^2))
    }
    
    x <- seq(from=1, to=2, length.out = 300)
    y <- lorentzian(x, 1.3, 0.01, 1) + lorentzian(x, 1.6, 0.01, 1)
    dataset <- new_nmr_dataset(
        metadata = list(external = data.frame(NMRExperiment = "10")),
        data_fields = list(
            data_1r = list(y)
        ),
        axis = list(list(x))
    )
    expect_warning(
        nmr_autophase(dataset, method="NLS"),
        "all_components=TRUE"
    )
})

test_that("nmr_autophase forwards imaginary-data availability to NMRphasing's absorptionOnly argument", {
    skip_if_not_installed("NMRphasing")

    # Force serial execution so the mocked binding (only visible in this
    # process) is actually used, and so results are recorded in a
    # deterministic, per-sample order.
    old_bpparam <- BiocParallel::bpparam()
    BiocParallel::register(old_bpparam)
    BiocParallel::register(BiocParallel::SerialParam())
    on.exit(BiocParallel::register(old_bpparam), add = TRUE)

    log_file <- tempfile()
    on.exit(unlink(log_file), add = TRUE)
    testthat::local_mocked_bindings(
        NMRphasing = function(X, absorptionOnly, ...) {
            # Record what nmr_autophase() actually asked for.
            cat(absorptionOnly, "\n", file = log_file, append = TRUE)
            X
        },
        .package = "NMRphasing"
    )

    x <- seq(from = 1, to = 2, length.out = 20)
    y_real <- x
    y_imag <- x * 2
    dataset <- new_nmr_dataset(
        metadata = list(external = data.frame(NMRExperiment = c("10", "20"))),
        data_fields = list(
            # Sample "10" has imaginary data, sample "20" does not.
            data_1r = list(y_real, y_real),
            data_1i = list(y_imag, NULL)
        ),
        axis = list(list(x), list(x))
    )

    expect_warning(
        nmr_autophase(dataset, method = "NLS"),
        "missing imaginary spectrum"
    )

    logged <- scan(log_file, what = character(), quiet = TRUE)
    expect_length(logged, 2)
    # Sample with imaginary data => absorptionOnly must be FALSE (full
    # complex phasing). Sample without it => absorptionOnly must be TRUE.
    expect_equal(logged, c("FALSE", "TRUE"))
})

test_that("nmr_autophase summarizes many missing-imaginary samples instead of listing them all", {
    skip_if_not_installed("NMRphasing")

    # NMRphasing::NMRphasing() is mocked here (as an identity function) purely
    # to keep the test fast and independent of the phasing algorithm's
    # numerical behavior on degenerate synthetic data -- we only care about
    # the message that is built before NMRphasing is invoked.
    testthat::local_mocked_bindings(
        NMRphasing = function(X, absorptionOnly, ...) X,
        .package = "NMRphasing"
    )

    old_width <- options(width = 1000)
    on.exit(options(old_width), add = TRUE)

    n <- 10
    x <- seq(from = 1, to = 2, length.out = 20)
    y <- x
    sample_names <- paste0("Sample", seq_len(n))
    dataset <- new_nmr_dataset(
        metadata = list(external = data.frame(NMRExperiment = sample_names)),
        data_fields = list(
            data_1r = stats::setNames(rep(list(y), n), sample_names),
            # None of the 10 samples has imaginary data.
            data_1i = stats::setNames(vector("list", n), sample_names)
        ),
        axis = rep(list(list(x)), n)
    )

    warnings_raised <- testthat::capture_warnings(
        nmr_autophase(dataset, method = "NLS")
    )
    combined <- paste(warnings_raised, collapse = " | ")

    # The truncated form must be used ...
    expect_match(combined, "and 5 more", fixed = TRUE)
    # ... and not the full dump of all 10 sample names (in particular the
    # 6th-to-10th names must not appear).
    for (missing_name in sample_names[6:10]) {
        expect_false(grepl(missing_name, combined, fixed = TRUE))
    }
})
