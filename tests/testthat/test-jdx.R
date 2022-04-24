test_that("comments are stripped", {
  
  expect_equal(strip_comments(c("##TITLE=Hello $$ world",
                                "$$ Hello world",
                                "12343 $$ hello world")),
               c("##TITLE=Hello ", "", "12343 "))
  
})

test_that("comments are stripped", {
  
  expect_equal(strip_comments(c("##TITLE=Hello $$ world",
                                "$$ Hello world",
                                "12343 $$ hello world")),
               c("##TITLE=Hello ", "", "12343 "))
  
})
