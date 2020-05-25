context("test-random_subsampling")

test_that("random_subsampling", {
subject_id <- c("Alice", "Bob", "Alice", "Bob")
replicate <- c(1, 1, 2, 2)
rnd <- random_subsampling(1:4, iterations = 2, test_size = 0.25, keep_together = subject_id)
rnd2 <- random_subsampling(1:4, iterations = 2, test_size = 0.25, keep_together = NULL)
rnd3 <- random_subsampling(1:3, iterations = 2, test_size = 0.3)
expect_true(is.list(rnd))
expect_true(is.list(rnd2))
expect_true(is.list(rnd3))
})
