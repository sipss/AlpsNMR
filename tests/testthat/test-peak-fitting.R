test_that("get_peak_bounds works", {
    # on the right, we are just between 3 and 4
    # on the left,  we are 80% closer to 2 than to 1
    expect_equal(
        get_peak_bounds(peak_limit_left = 1, peak_limit_right = 3, pos = 2, x = 1:4, sgf = c(-4, 1, 2, -2)),
        c(left = 1, apex = 2, right = 3, xleft = 1.8, xright = 3.5)
    )
})
