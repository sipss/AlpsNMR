# Regression tests for bugs found and fixed in a prior critical code review
# of R/area_estimation.R. Each test below targets one specific fix.


# ---------------------------------------------------------------------------
# 1. get_norm_rmse() argument-collision fix
# ---------------------------------------------------------------------------
# get_norm_rmse(y_fitted, y, y_fitted_apex, y_apex) used to be called at its
# only call site (inside peaklist_fit_lorentzians()) with a mix of positional
# and named arguments that collided by exact name-tag matching, silently
# scrambling which vector was bound to which parameter. The call site is now
# fixed to pass all 4 arguments explicitly by name (see the call around
# `norm_rmse <- get_norm_rmse(...)` in R/area_estimation.R).
#
# get_norm_rmse() itself is a small pure function of its 4 arguments, so we
# can directly verify its formula/contract here:
#   y_fitted_adj <- y_fitted - y_fitted_apex + y_apex
#   mse          <- sum((y - y_fitted_adj)^2) / length(y)
#   rmse         <- sqrt(mse)
#   norm_rmse    <- rmse / y_fitted_apex
# This test does not (and cannot) directly re-detect a positional-argument
# regression at the call site -- that is instead guarded by confirming below
# that the call site source text uses all-named arguments.
test_that("get_norm_rmse() computes the documented formula from explicitly-named arguments", {
    y_fitted <- c(1, 2, 3, 4, 5)
    y <- c(1.1, 2.1, 2.9, 4.2, 4.8)
    y_fitted_apex <- 3
    y_apex <- 2.9

    result <- get_norm_rmse(
        y_fitted = y_fitted,
        y = y,
        y_fitted_apex = y_fitted_apex,
        y_apex = y_apex
    )

    y_fitted_adj <- y_fitted - y_fitted_apex + y_apex
    expected_mse <- sum((y - y_fitted_adj)^2) / length(y)
    expected_rmse <- sqrt(expected_mse)
    expected_norm_rmse <- expected_rmse / y_fitted_apex

    expect_equal(result, expected_norm_rmse)
    # Sanity: for these inputs the adjustment and rmse are non-trivial
    # (i.e. this isn't accidentally testing a degenerate 0 == 0 case).
    expect_true(expected_rmse > 0)
})

test_that("get_norm_rmse() output changes when the apex arguments are swapped", {
    # If y_fitted_apex and y_apex were ever accidentally swapped (the exact
    # failure mode of the original bug), the result must differ -- this
    # confirms the function is actually sensitive to which value is which,
    # so a future regression at the call site would (in principle) be
    # detectable by comparing against the correct formula.
    y_fitted <- c(2, 4, 6, 8, 10)
    y <- c(2.5, 4.2, 5.7, 8.3, 9.6)

    result_correct <- get_norm_rmse(y_fitted = y_fitted, y = y, y_fitted_apex = 6, y_apex = 5.7)
    result_swapped <- get_norm_rmse(y_fitted = y_fitted, y = y, y_fitted_apex = 5.7, y_apex = 6)

    expect_false(isTRUE(all.equal(result_correct, result_swapped)))
})

test_that("peaklist_fit_lorentzians() calls get_norm_rmse() with all 4 arguments named", {
    # Direct guard against the exact regression: re-introducing a positional
    # (or partially-named) call at the call site. We inspect the actual
    # function body rather than the file text, so this survives reformatting.
    fn_body <- deparse(body(peaklist_fit_lorentzians))
    call_line_idx <- grep("get_norm_rmse\\(", fn_body)
    expect_true(length(call_line_idx) == 1)

    # Grab the call expression from the parsed function body directly.
    find_call <- function(expr) {
        if (is.call(expr)) {
            if (identical(as.character(expr[[1]]), "get_norm_rmse")) {
                return(expr)
            }
            for (part in as.list(expr)[-1]) {
                # Some sub-expressions (e.g. empty/missing arguments in
                # indexing calls) raise "argument is missing" merely by
                # being referenced; skip those defensively.
                found <- tryCatch(find_call(part), error = function(e) NULL)
                if (!is.null(found)) {
                    return(found)
                }
            }
        }
        NULL
    }
    call_expr <- find_call(body(peaklist_fit_lorentzians))
    expect_false(is.null(call_expr))

    arg_names <- names(as.list(call_expr))[-1] # drop the function name itself
    expect_equal(sort(arg_names), sort(c("y_fitted", "y", "y_fitted_apex", "y_apex")))
    expect_true(all(nzchar(arg_names)))
})


# ---------------------------------------------------------------------------
# 2. get_peak_bounds() NA sentinel fix
# ---------------------------------------------------------------------------
# get_peak_bounds(peak_limit_left, peak_limit_right, pos, x, sgf) used to
# return c(left=0, right=0, xleft=0, xright=0) when no inflection point
# existed on one side of the peak (e.g. peak too close to the spectrum
# edge), a 0 sentinel indistinguishable from a genuine ppm/index value.
# It must now return NA_real_ for left/right/xleft/xright in that case.

test_that("get_peak_bounds() returns NA (not 0) when there is no left inflection point", {
    # peak_limit_left = 5 with pos = 2 is filtered out by
    # `peak_limit_left[pos - peak_limit_left > 0]` (2 - 5 = -3, not > 0),
    # leaving zero candidate left inflection points -- the edge-of-spectrum
    # case for a peak near the start of the spectrum.
    result <- get_peak_bounds(
        peak_limit_left = 5,
        peak_limit_right = 8,
        pos = 2,
        x = 1:10,
        sgf = rep(0, 10)
    )

    expect_true(is.na(result[["left"]]))
    expect_true(is.na(result[["right"]]))
    expect_true(is.na(result[["xleft"]]))
    expect_true(is.na(result[["xright"]]))
    expect_equal(result[["apex"]], 2)

    # It must be a real NA_real_, not a numeric 0 masquerading as missing.
    expect_false(isTRUE(result[["left"]] == 0))
})

test_that("get_peak_bounds() returns NA (not 0) when there is no right inflection point", {
    # peak_limit_right is empty outright -- the edge-of-spectrum case for a
    # peak near the end of the spectrum.
    result <- get_peak_bounds(
        peak_limit_left = 1,
        peak_limit_right = integer(0),
        pos = 9,
        x = 1:10,
        sgf = rep(0, 10)
    )

    expect_true(is.na(result[["left"]]))
    expect_true(is.na(result[["right"]]))
    expect_true(is.na(result[["xleft"]]))
    expect_true(is.na(result[["xright"]]))
    expect_equal(result[["apex"]], 9)
})

test_that("get_peak_bounds() still returns real (non-NA) bounds in the normal case", {
    # Same fixture as the existing get_peak_bounds() test in
    # test-peak-fitting.R, confirming the fix didn't disturb the happy path.
    result <- get_peak_bounds(
        peak_limit_left = 1,
        peak_limit_right = 3,
        pos = 2,
        x = 1:4,
        sgf = c(-4, 1, 2, -2)
    )

    expect_equal(
        result,
        c(left = 1, apex = 2, right = 3, xleft = 1.8, xright = 3.5)
    )
    expect_false(anyNA(result))
})


# ---------------------------------------------------------------------------
# 3. refine_lorentzian_fit_with_nls() error-capture fix
# ---------------------------------------------------------------------------
# The tryCatch() error= handler used to assign to a LOCAL `error_msgs`
# inside its own closure (missing the `<<-` superassignment), so nls
# convergence/evaluation errors were silently swallowed and
# new_params[["error_msgs"]] was always NULL. It now uses `<<-`, so the
# error message from a failing stats::nls() call is captured and returned.

test_that("refine_lorentzian_fit_with_nls() captures nls errors via error_msgs", {
    # Two data points cannot identify a 3-parameter lorentzian model
    # (A, x0, gamma): stats::nls() fails immediately (e.g. "missing value
    # or an infinity produced when evaluating the model" / singular
    # gradient), which is exactly the kind of error this function must
    # capture and report rather than discard.
    data_to_fit <- data.frame(x = c(1, 2), y = c(1, 2))
    start <- list(A = 1, x0 = 1.5, gamma = 0.1)

    result <- refine_lorentzian_fit_with_nls(data_to_fit, start, method = "peak")

    expect_false(is.null(result[["error_msgs"]]))
    expect_true(is.character(result[["error_msgs"]]))
    expect_true(nchar(result[["error_msgs"]][1]) > 0)

    # When nls fails, the function must fall back to the starting values
    # rather than propagating the error to the caller.
    expect_equal(result[["estimated_A"]], start[["A"]])
    expect_equal(result[["gamma"]], start[["gamma"]])
    expect_equal(result[["peak_pos_ppm"]], start[["x0"]])
})

test_that("refine_lorentzian_fit_with_nls() leaves error_msgs empty when nls succeeds", {
    # Sanity/contrast case: a clean, well-sampled lorentzian (plus a tiny
    # amount of noise, since a perfectly noiseless exact fit can itself
    # confuse nls's relative-offset convergence criterion) should let nls
    # converge, so no error should be recorded.
    set.seed(42)
    lorentzian <- function(x, x0, gamma, A) {
        A * (1 / (pi * gamma)) * ((gamma^2) / ((x - x0)^2 + gamma^2))
    }
    x <- seq(-5, 5, length.out = 200)
    true_x0 <- 0.3
    true_gamma <- 1.1
    true_A <- 4
    y <- lorentzian(x, true_x0, true_gamma, true_A) + rnorm(length(x), sd = 0.0005)
    data_to_fit <- data.frame(x = x, y = y)
    start <- list(A = true_A * 0.9, x0 = 0, gamma = true_gamma * 0.9)

    result <- refine_lorentzian_fit_with_nls(data_to_fit, start, method = "peak")

    expect_null(result[["error_msgs"]])
    expect_equal(result[["peak_pos_ppm"]][[1]], true_x0, tolerance = 1e-2)
})


# ---------------------------------------------------------------------------
# 4. peaklist_fit_lorentzians() all_errors$error_msg accumulation fix
# ---------------------------------------------------------------------------
# Previously `all_errors$error_msg <- paste(new_params[["error_msgs"]],
# collapse = "\n")` OVERWROTE the error_msg vector on every loop iteration,
# while `all_errors$peak_id` grew via c() -- so after 2+ peaks hit an error
# path, `peak_id` and `error_msg` ended up different lengths (and
# as.data.frame() on a list with mismatched-length elements either warns/
# recycles or errors). It's now grown via c() in lockstep with peak_id.
#
# To trigger this without a full nmr_dataset pipeline, we build a minimal
# 1-sample nmr_dataset_1D with two well-separated, clean lorentzian peaks
# (so get_peak_bounds() finds real, non-NA bounds via the raw signal's
# inflection points) but supply a jagged, non-lorentzian "baseline"
# (data_1r_baseline) inside each peak's fitting window. Since
# refine_peak_model = "peak" fits the *baseline-corrected* data, the nls
# fit is handed non-lorentzian-shaped data for both peaks and fails to
# converge for both -- forcing both peaks down the error-accumulating
# branch.
test_that("peaklist_fit_lorentzians() accumulates peak_id and error_msg in lockstep", {
    lorentzian <- function(x, x0, gamma, A) {
        A * (1 / (pi * gamma)) * ((gamma^2) / ((x - x0)^2 + gamma^2))
    }

    npoints <- 80
    x <- seq(0, 10, length.out = npoints)
    dx <- x[2] - x[1]

    y <- rep(0.01, npoints) +
        lorentzian(x, x[20], gamma = dx * 3, A = 5) +
        lorentzian(x, x[60], gamma = dx * 3, A = 5)

    # A jagged "baseline" confined to each peak's fitting window (13:26 and
    # 53:66, as determined by the inflection points of `y` above) so that
    # y - y_basel is not lorentzian-shaped there and nls fails to converge.
    y_basel <- rep(0, npoints)
    jag <- function(idxs) rep(c(50, -50), length.out = length(idxs))
    y_basel[13:26] <- jag(13:26)
    y_basel[53:66] <- jag(53:66)

    nmr_dataset <- new_nmr_dataset_1D(
        ppm_axis = x,
        data_1r = matrix(y, nrow = 1),
        metadata = list(external = data.frame(NMRExperiment = "10"))
    )
    nmr_dataset$data_1r_baseline <- matrix(y_basel, nrow = 1)

    peak_data <- data.frame(
        peak_id = c("Peak1", "Peak2"),
        NMRExperiment = c("10", "10"),
        ppm = c(x[20], x[60]),
        pos = c(20, 60),
        intensity = c(y[20], y[60])
    )

    withCallingHandlers(
        {
            result <- peaklist_fit_lorentzians(
                peak_data,
                nmr_dataset,
                amplitude_method = "intensity",
                refine_peak_model = "peak"
            )
        },
        warning = function(w) {
            if (grepl("longer object length", conditionMessage(w))) {
                fail(paste("Recycling warning (length mismatch bug reappeared):", conditionMessage(w)))
            }
            invokeRestart("muffleWarning")
        }
    )

    errors <- attr(result, "errors")
    expect_true(is.data.frame(errors))
    expect_true(nrow(errors) >= 2)
    expect_equal(length(errors$peak_id), length(errors$error_msg))
    expect_setequal(errors$peak_id, c("Peak1", "Peak2"))
    expect_true(all(nzchar(errors$error_msg)))
})
