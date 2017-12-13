context("Import JDX files")

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

context("list_of_lists_to_tibble")

test_that("list_of_lists_to_tibble works", {
  data_as_list <- list(list(Gender = "Male", Height = 170),
                      list(Gender = "Female", Height = 160),
                      list(Gender = "Male", Weight = 80))
  mydata <- list_of_lists_to_tibble(data_as_list)
  # check result is consistent:
  expect_equal(nrow(mydata), 3)
  expect_equal(mydata[["Gender"]], c("Male", "Female", "Male"))
  expect_equal(mydata[["Height"]], c(170, 160, NA))
  expect_equal(mydata[["Weight"]], c(NA, NA, 80))

})

test_that("list_of_lists_to_tibble keeps NULL rows as NA rows", {
  data_as_list <- list(list(Gender = "Male", Height = 170),
                       list(Gender = "Female", Height = 160),
                       NULL,
                       list(Gender = "Male", Weight = 80))
  mydata <- list_of_lists_to_tibble(data_as_list)
  # check result is consistent:
  expect_equal(nrow(mydata), 4)
  expect_equal(mydata[["Gender"]], c("Male", "Female", NA, "Male"))
  expect_equal(mydata[["Height"]], c(170, 160, NA, NA))
  expect_equal(mydata[["Weight"]], c(NA, NA, NA, 80))
})

test_that("list_of_lists_to_tibble replaces NULL values with NA values", {
  data_as_list <- list(list(Gender = "Male", Height = 170),
                       list(Gender = "Female", Height = NULL),
                       NULL,
                       list(Gender = "Male", Weight = 80))
  mydata <- list_of_lists_to_tibble(data_as_list)
  # check result is consistent:
  expect_equal(nrow(mydata), 4)
  expect_equal(mydata[["Gender"]], c("Male", "Female", NA, "Male"))
  expect_equal(mydata[["Height"]], c(170, NA, NA, NA))
  expect_equal(mydata[["Weight"]], c(NA, NA, NA, 80))
})

test_that("list_of_lists_to_tibble accepts POSIXct", {
  data_as_list <- list(list(Gender = "Male", Height = 170, Day = as.POSIXct("2010-01-01")),
                       list(Gender = "Female", Height = NULL, Day = NULL),
                       NULL,
                       list(Gender = "Male", Weight = 80, Day = NA))
  mydata <- list_of_lists_to_tibble(data_as_list)
  # check result is consistent:
  expect_equal(nrow(mydata), 4)
  expect_equal(mydata[["Gender"]], c("Male", "Female", NA, "Male"))
  expect_equal(mydata[["Height"]], c(170, NA, NA, NA))
  expect_equal(mydata[["Weight"]], c(NA, NA, NA, 80))
  expect_equal(mydata[["Day"]], as.POSIXct(c("2010-01-01", NA, NA, NA)))
})


context("simplify_df simplifies a data frame")

test_that("simplify_df works", {
  mydata <- data.frame(interesting = 1:5,
                       stupid = rep(1, 5),
                       interesting2 = 2*(1:5))
  output <- simplify_df(mydata)
  expect_equal(class(output), "list")
  expect_equal(names(output), c("const", "diff"))
  expect_equal(output[["const"]], data.frame(stupid = 1))
  expect_equal(output[["diff"]], mydata[,c("interesting", "interesting2"), drop = FALSE])
})

