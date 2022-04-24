test_that("random_subsampling", {
subject_id <- c("Alice", "Bob", "Chalie", "Eve")
rnd <- random_subsampling(1:4, iterations = 2, test_size = 0.25, keep_together = subject_id)
rnd2 <- random_subsampling(1:4, iterations = 2, test_size = 0.25, keep_together = NULL)
rnd3 <- random_subsampling(1:4, iterations = 2, test_size = 0.3)
expect_true(is.list(rnd))
expect_true(is.list(rnd2))
expect_true(is.list(rnd3))
})
